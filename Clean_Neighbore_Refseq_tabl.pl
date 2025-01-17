#!/usr/bin/env perl
use strict;
use warnings;
use Modern::Perl '2011';
use autodie;
use Smart::Comments;

# Validate arguments
unless (@ARGV == 1) {
    die <<"EOT";
Usage: $0 <infile>
This software cleans virus neighbor data for non-human viruses.
Example: $0 Human_Viruses_RefSeq_and_neighbors_genome_data.tab
EOT
}

# Input and output files
my $infile = shift;
my ($outname) = $infile =~ m/(.*)\.\w+/;
my $outfile = "$outname-cleaned.tsv";

open my $in, '<', $infile;
open my $out, '>', $outfile;

# Exclusion patterns
my @patterns = qw(
    delavirus phage satellite iridovirus entomopox myoviridae monkey
    lemur gorilla chimpanzee swine powassan chicken marseillevirus
    ilarvirus kappatorquevirus alphacarmotetravirus orthobunyavirus
    iotatorquevirus potexvirus giardiavirus macavirus baculovirus
    gammacarmovirus macaque ranid cercopithecine canary taterapox
    fowlpox macaca mimiviridae feline siphoviridae frog tanapox
    phycodnaviridae tevenvirinae betabaculovirus cyprinid buffalo
    ectromelia molluscum bornavirus elephant salmon simian mammarenavirus
    gammaretrovirus cuevavirus henipavirus orthohantavirus marburgvirus
    baculoviridae porcine equine mycobacterium cafeteriavirus bat
    crocodile muromegalo toursvirus betasatellite ranavirus mollusc
    probosci versivirus avian reptarenavirus aviadeno d3112 leishmania
    g4micr l5virus iflavi potyvirus ipomovirus caulimo murine spuma
    atadeno likavirus suipox autographi canine t4virus avula lago
    bornavirus centapox chordopo clostero alphabaculo kp282679 grablo
    protoparv coccolit chlorovirus microviridae asfarvi bovine scuta
    begomo avastro betaentomo kc977571 podoviri sripu badna phlebo
    marafi cyprini tymovirus camel horse murray pteropin foxe whale
    pigeon alouatta alpaca btpa-betacov ateles btvs-betacov chandipura
    coronavirus neoromicia crimean-congo cyclovirus duck echovirus eilat
    equid ferret fitzroy gb gemycircularvirus jurona kamiti kemerovo
    lordsdale macacine mammalian onyong-nyong panine phopivirus
    piliocolobus
);

# Combine patterns into a single regex
my $exclude_regex = join '|', map { quotemeta } @patterns;

# Process file
while (my $line = <$in>) {
    chomp $line;
    my $lc_line = lc $line;

    # Skip lines matching exclusion patterns
    next if $lc_line =~ /$exclude_regex/xms;

    # Write remaining lines to output
    say $out $line;
}

close $in;
close $out;
