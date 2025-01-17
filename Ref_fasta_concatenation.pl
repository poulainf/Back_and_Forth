#!/usr/bin/env perl
use warnings;
use strict;
use Modern::Perl '2011';
use autodie;
use Smart::Comments;
use warnings FATAL => 'uninitialized';
use List::MoreUtils qw(uniq);
use List::Util qw/shuffle/;

unless (@ARGV == 3) {
    die << "EOT";
Usage: $0 <infile.fasta> <infile.gb> <infile.tsv>
This tool reformats FASTA IDs for TimeTree and BEAST analysis.
Example: $0 test.fasta test.gb test.tsv
EOT
}

### Input files
my $seq_file  = shift;
my $info_file = shift;
my $list_file = shift;

my ($outname) = $seq_file =~ m/(.*)\.\w+/;

### Sequence analysis
my %seq_for = read_fasta($seq_file);
my %list_for = list_id($list_file);

## Process information
my ($info_for_ref, $info_carac_ref) = taxon_id($info_file, %list_for);
my %info_for   = %$info_for_ref;
my %info_carac = %$info_carac_ref;

my @criteria = qw(tax_file taxon genotype subtype serotype group type lineage);

### Start processing
for my $species (keys %info_for) { ## Elapsed time |===[%]
    for my $carac (@criteria) {
        for my $taxa (keys %{ $info_carac{$species}{$carac} }) {
            next if $taxa eq "NA";

            my @IDs = @{ $info_carac{$species}{$carac}{$taxa} };
            if (@IDs > 9) {
                my $taxo = $taxa;
                $taxo =~ s{/}{-}g;

                my $prefix = "NO-iT";
                $prefix = "OK-iT" if @IDs > 49;

                my $outfile = join "_", $prefix, $species, $carac, "$taxo.fa";
                open my $out, '>', $outfile;

                for my $ID (@IDs) {
                    my $date = $info_for{$species}{$ID};
                    if ($seq_for{$ID}) {
                        say {$out} join "|", ">$ID", $date;
                        say {$out} $seq_for{$ID};
                    }
                }
                close $out;
            }
        }
    }
}

### Functions

sub read_fasta {
    my $infile = shift;
    open my $in, '<', $infile;

    my (%seq_for, $seq_id, $seq);

    while (my $line = <$in>) {
        chomp $line;
        if ($line =~ m/^>([^.]+)\..*/) {
            $seq_for{$seq_id} = $seq if $seq_id;
            $seq_id = $1;
            $seq = '';
        } else {
            $seq .= $line;
        }
    }
    $seq_for{$seq_id} = $seq if $seq_id;
    close $in;

    return %seq_for;
}

sub list_id {
    my $infile = shift;
    open my $in, '<', $infile;

    my %list_for;

    while (my $line = <$in>) {
        chomp $line;
        next if $line =~ m/^#/;

        my @fields = split /\t/, $line;
        my ($ID1, $ID, $host_name, $specie_name, $seg) = @fields[0, 1, 2, 4, 5];
        $host_name = lc $host_name;
        $specie_name = lc $specie_name;
        $seg = lc $seg;

        $specie_name =~ s/human//;
        $specie_name =~ s/\s+$//;
        $specie_name =~ s/^\s+//;
        $specie_name =~ tr/ /_/;

        $seg =~ s/\s+$//;
        $seg =~ s/^\s+//;
        $seg =~ tr/ /_/;

        my $tax = join "_", $specie_name, $seg;

        $list_for{$ID} = {
            host_name => $host_name,
            seg       => $seg,
            tax       => $tax
        };

        for my $sub_id (split /,/, $ID1) {
            $list_for{$sub_id} = {
                host_name => $host_name,
                seg       => $seg,
                tax       => $tax
            };
        }
    }
    close $in;
    return %list_for;
}

sub taxon_id {
    my ($infile, %list_for) = @_;
    open my $in, '<', $infile;

    my %month_ref = (
        jan => "01", feb => "02", mar => "03", apr => "04",
        may => "05", jun => "06", jul => "07", aug => "08",
        sep => "09", oct => "10", nov => "11", dec => "12",
    );

    my (%info_for, %info_carac);
    my ($ID, $species, $time, $taxon, $host, $genotype, $subtype, $serotype, $group, $type, $lineage);
    my $next_line = "OFF";

    while (my $line = <$in>) {
        chomp $line;
        next if $line =~ m/mol_type=|rpt_type|TITLE|Phenotype/x;

        if ($line =~ m/VERSION\s+([^.]+)/) {
            ($ID, $time, $taxon, $host, $genotype, $subtype, $serotype, $group, $type, $lineage) = ($1, "None", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA");
            $next_line = "OFF";
        } elsif ($line =~ m/ORGANISM\s+(.+)/) {
            $species = lc $1;
        } elsif ($line =~ m/collection_date=['"]?(.+)["']?/) {
            $time = $1;
        } elsif ($line =~ m/^\/\/$/) {
            # Process collected data
            $time = format_date($time, \%month_ref);
            my $tax = $list_for{$ID}{tax};
            my $seg = $list_for{$ID}{seg};
            next unless $time ne "None" && $host =~ m/human|homo/i;

            $info_for{$tax}{$ID} = $time;
            push @{ $info_carac{$tax}{taxon}{$taxon} }, $ID if $taxon ne "NA";
            # Add additional fields to info_carac as needed
        }
    }
    close $in;
    return (\%info_for, \%info_carac);
}

sub format_date {
    my ($time, $month_ref) = @_;
    # Add logic to handle and standardize date formats
    return $time;
}
