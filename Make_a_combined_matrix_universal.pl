#!/usr/bin/env perl
use warnings;
use strict;
use Modern::Perl '2011';
use autodie;
#use Smart::Comments;
use Algorithm::Numerical::Shuffle qw/shuffle/;
use List::Util qw(sum);
use List::Util qw(max);
use List::UtilsBy qw(max_by);
use warnings FATAL => 'uninitialized';
use Tie::IxHash;


### Input files
my $file = shift;
my $the_in = shift;
my $rounding =shift;
my $virus;


if ($file =~ m/.*(\w\w\-iT.+N\d+.+)_ALIGNED.+matrix.*\.txt/xms){ 
	$virus=$1;
	}
else { $virus=$1 if $file =~ m/(.+)\.txt/xms};


my %mix_matrix;
tie %mix_matrix, 'Tie::IxHash';

if ($rounding eq "Rate" ){
	
	%mix_matrix = read_fasta_rate( $file );

	}else{

	%mix_matrix = read_fasta( $file );

};

my @values_hash= keys %mix_matrix;
my $first_key = $values_hash[0];






## %mix_matrix

my $outfile_cons;
if ($the_in){
	$outfile_cons = join "_" ,"Combined",$the_in,$virus,'mediane_matrix.txt';
}
else{
	$outfile_cons = join "_" ,"Combined",$virus,'mediane_matrix.txt';
};




	
open my $ouit, '>',  $outfile_cons;
say {$ouit} join "\t", "Mutation\ Types",$virus;

for my $tres(keys %mix_matrix){
	
	say {$ouit} join "\t", $tres, @{$mix_matrix{$tres}};

};



# Functions

sub median{
	
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}


sub read_fasta {

	my %mix_matrix;
	tie %mix_matrix, 'Tie::IxHash';
	
	my $file = shift ;
			
	open my $in, '<', $file;
		## $virus
	LINE:
	while (my $line = <$in> ){
		
		chomp $line;
		if ($line =~ m/Mutation\ Types/xms){next LINE};
		if ($line =~ m/\#/xms){next LINE};
		## $line
		my @values = split "\t" ,$line;
		my $remove = shift @values;				
		my $mediane= median(@values);
		my $rounded_mediane = sprintf "%.0f", $mediane;
		push @{$mix_matrix{$remove}}, $rounded_mediane;
	
	}
		
	close $in;
	return %mix_matrix;	

}



sub read_fasta_rate {

	my %mix_matrix;
	tie %mix_matrix, 'Tie::IxHash';
	
	my $file = shift ;
			
	open my $in, '<', $file;
		## $virus
	LINE:
	while (my $line = <$in> ){
		
		chomp $line;
		if ($line =~ m/Mutation\ Types/xms){next LINE};
		if ($line =~ m/\#/xms){next LINE};
		## $line
		my @values = split "\t" ,$line;
		my $remove = shift @values;
				
				
				
		my $mediane= median(@values);
		#my $rounded_mediane = sprintf "%.0f", $mediane;

		#push @{$mix_matrix{$remove}}, $rounded_mediane;
		push @{$mix_matrix{$remove}}, $mediane;
	
	}
		
	close $in;
	return %mix_matrix;	

}
