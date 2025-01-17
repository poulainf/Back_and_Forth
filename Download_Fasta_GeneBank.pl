#!/usr/bin/env perl
use warnings;
use strict;
use Modern::Perl '2011';
use autodie;
use Smart::Comments;
use warnings FATAL => 'uninitialized';

### Input files
my $seq_file = shift or die "Usage: $0 <sequence_file>\n";
my %Datas = read_fasta($seq_file);

TOUR:
for my $id (keys %Datas) { ### Elapsed time |===[%]

    my $specie_name = $Datas{$id}{"Taxo"};
    my $segment     = $Datas{$id}{"Segment"};
    my @list_ids    = @{$Datas{$id}{"List"}};

    # Skip if no IDs available
    if ($list_ids[0] eq 'NULL') {
        next TOUR;
    }

    # Include the current ID in the list
    push @list_ids, $id;

    # Create a directory for the ID if it doesn't exist
    my $comm = "[ ! -d $id ] && mkdir -p $id";
    system($comm) == 0 or warn "Failed to create directory for $id: $!\n";

    TORN:
    while (@list_ids) {
        my $name = shift @list_ids;

        # Construct filename
        my $filename = "./$id/$name.fast";

        # Fetch data only if file doesn't exist or is empty
        unless (-e $filename && -s $filename) {
            my $command = "esearch -db Nucleotide -query \"$name\" < /dev/null | efetch -format fasta > \"$filename\"";
            system($command) == 0 or warn "Failed to fetch FASTA for $name: $!\n";
            next TORN;
        }

    }

    next TOUR;
}

### Subroutine to read FASTA data
sub read_fasta {
    my $infile = shift or die "Input file not specified for read_fasta\n";

    ## Reading input file: $infile
    ## Elapsed time |===[%]

    open my $in, '<', $infile;

    my %dna;

    TURN:
    while (my $line = <$in>) {
        chomp $line;

        # Skip comment lines
        if ($line =~ m/^#/xms) {
            next TURN;
        }

        # Split the line into fields
        my @line = split "\t", $line;
        ## @line

        my ($family, $genus, $specie);

        # Extract fields from line
        my $Neighbor       = $line[1];
        my $Representative = shift @line;
        my $Segment        = pop @line;
        my $Taxonomy       = pop @line;
        my $Selected       = pop @line;
        my $Host           = pop @line;

        # Parse taxonomy field
        if ($Selected =~ m/([^,]+),([^,]+),([^,]+)/xms) {
            $family = $1;
            $genus  = $2;
            $specie = $3;
        }

        # Populate DNA hash
        if (exists $dna{$Representative}) {
            push @{$dna{$Representative}{"List"}}, $Neighbor;
        } else {
            $dna{$Representative} = {
                List    => \@line,
                Segment => $Segment,
                Taxo    => $Taxonomy,
                Lineage => $Selected,
                Host    => $Host,
            };
        }

        next TURN;
    }

    close $in;
    return %dna;
}
