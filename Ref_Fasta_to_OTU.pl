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
This tool reformats FASTA IDs for TreeTime and BEAST analysis.
Example: $0 test.fasta test.gb test.tsv
EOT
}

### Input files
my $seq_file  = shift;
my $info_file = shift;
my $list_file = shift;

my ($outname) = $seq_file =~ m/(.*)\.\w+/;

### Sequence analysis
my %seq_for   = read_fasta($seq_file);
my %list_for  = list_id($list_file);
my ($info_for_ref, $info_carac_ref) = taxon_id($info_file, %list_for);
my %info_for  = %$info_for_ref;
my %info_carac = %$info_carac_ref;

my @criteria = ("tax_file", "taxon", "genotype", "subtype", "serotype", "group", "type", "lineage");

### Start of processing
for my $species (keys %info_for) {    ## Elapsed time |===[%]
    for my $carac (@criteria) {
        for my $taxa (keys %{ $info_carac{$species}{$carac} }) {
            next if $taxa eq "NA";

            my @IDs = @{ $info_carac{$species}{$carac}{$taxa} };

            if (scalar @IDs > 9) {
                my $taxo = $taxa;
                $taxo =~ s/\//-/g;

                my $prefix = "NO-iT";
                $prefix = "OK-iT" if scalar @IDs > 49;

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

    my %seq_for;
    my ($seq_id, $seq);

    while (my $line = <$in>) {
        chomp $line;
        if ($line =~ m/^>([^\.]+)\..*/) {
            $seq_for{$seq_id} = $seq if $seq_id;
            $seq_id = $1;
            $seq    = '';
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

        $host_name  = lc $host_name;
        $specie_name = lc($specie_name);
        $specie_name =~ s/human//g;
        $specie_name =~ s/^\s+|\s+$//g;
        $specie_name =~ tr/ /_/;

        $seg =~ s/^\s+|\s+$//g;
        $seg =~ tr/ /_/;

        my $tax = join "_", $specie_name, $seg;

        $list_for{$ID} = {
            host_name => $host_name,
            seg       => $seg,
            tax       => $tax
        };

        for my $ID2 (split /,/, $ID1) {
            $list_for{$ID2} = $list_for{$ID};
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
        sep => "09", oct => "10", nov => "11", dec => "12"
    );

    my (%info_for, %info_carac);

    my ($ID, $specie, $time, $taxon, $host, $genotype, $subtype, $serotype, $group, $type, $lineage);
    $time = $taxon = $host = $genotype = $subtype = $serotype = $group = $type = $lineage = "NA";

    while (my $line = <$in>) {
        chomp $line;

        if ($line =~ m/^VERSION\s+([^\s]+)/) {
            $ID      = $1;
            $time    = $taxon = $host = $genotype = $subtype = $serotype = $group = $type = $lineage = "NA";
        } elsif ($line =~ m/^ +ORGANISM +(.+)/) {
            $specie = lc $1;
        } elsif ($line =~ m/collection_date=["']?([^"']+)["']?/) {
            $time = $1;
        } elsif ($line =~ m/^\/\//) {
            $time =~ s#/#-#g;
            $time =~ s/^(\d{2})-(\w+)-(\d{4})$/$3-$month_ref{$2}-$1/;

            my $segment    = $list_for{$ID}{"seg"};
            my $tax        = $list_for{$ID}{"tax"};
            my $host_name  = $list_for{$ID}{"host_name"};

            next if exists $info_for{$tax}{$ID};
            $info_for{$tax}{$ID} = $time;

            push @{ $info_carac{$tax}{"taxon"}{$taxon} }, $ID;
            push @{ $info_carac{$tax}{"genotype"}{$genotype} }, $ID;
            push @{ $info_carac{$tax}{"subtype"}{$subtype} }, $ID;
            push @{ $info_carac{$tax}{"serotype"}{$serotype} }, $ID;
            push @{ $info_carac{$tax}{"group"}{$group} }, $ID;
            push @{ $info_carac{$tax}{"type"}{$type} }, $ID;
            push @{ $info_carac{$tax}{"lineage"}{$lineage} }, $ID;
        }
    }

    close $in;
    return (\%info_for, \%info_carac);
}
