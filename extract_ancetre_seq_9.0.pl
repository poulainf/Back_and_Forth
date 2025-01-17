#!/usr/bin/env perl
use warnings;
use strict;
use Modern::Perl '2011';
use autodie;
use Smart::Comments;
use List::AllUtils 'mesh';
use Algorithm::Numerical::Shuffle qw/shuffle/;
use List::Util qw(sum);
use List::Util qw(max);
use List::Util qw(min);
use warnings FATAL => 'uninitialized';
use Text::Balanced qw (
    extract_delimited
    extract_bracketed
    extract_quotelike
    extract_codeblock
    extract_variable
    extract_tagged
    extract_multiple
    gen_delimited_pat
    gen_extract_tagged
);
use Regexp::Common qw(balanced);
use Tie::IxHash;
#use Statistics::Basic qw(:all nofill);

my $tree_file = shift;
my $ID_file = shift;
my $INFO_table= shift;
my $INFO_size= shift;


my ($outname) = $tree_file =~ m/(.*)\.\w+/g;
my ($outname_short) = $tree_file =~ m/Uniq\_true\_tree\_(.+)\_sub_selected\_iteration\_\d+\.tree/g;

my ($ID_ref,$Id_ref_rev) = read_ID( $ID_file );
my %Id_ref_rev = %$Id_ref_rev;
my %ID_ref = %$ID_ref;

## %ID_ref
my %Info_hash = read_info( $INFO_table);
my ($tree_hash, $result_hash, $ID_feuilles) = read_tree( $tree_file);
my %tree_hash = %$tree_hash;
my %result_hash = %$result_hash;
my %ID_feuilles = %$ID_feuilles; 


my %Info_length = read_length( $INFO_size);


### $outname
### $outname_short

my $the_seq_length = $Info_length{$outname_short} ;

### $the_seq_length

my @tmp_length = keys %Info_length;

if (not defined $the_seq_length){

	my $tarttuf = $outname_short;
	$tarttuf =~ s/\_iter\-N\d+.+//g;
	my @tmp_grepage = grep (/${tarttuf}/, @tmp_length);
	my $tmp_grepage = $tmp_grepage[0];
	
	### $the_seq_length
	
	if (defined $tmp_grepage){
		
		### AAAA
		### $the_seq_length
	
		$the_seq_length = $Info_length{$tmp_grepage};
	
	}else{
		
	#	$tarttuf =~ s/\w\w\-iT\_//g;
	#	my @speciespecies = split "_", $tarttuf;
		
	#	my $tarttuf2 = "TATA";
		
	#	LENELENELEN:
	#	while (not defined $the_seq_length){
			
	#		$tarttuf =~ s/\_[^\_]*?$//;
			
			
	#		my @tmp_grepage_recherch = grep (/${tarttuf}/, @tmp_length);
			
	#		my $tmp_grepage_recherch = $tmp_grepage_recherch[0];

	#		if (defined $tmp_grepage_recherch){
	
	#			$the_seq_length = $Info_length{$tmp_grepage_recherch};
	
	#		};
			 
	#		if ($tarttuf eq $tarttuf2){ last LENELENELEN };
	#		$tarttuf2 = $tarttuf;
			
	#		next LENELENELEN;
		
	#	};
		
		
		if (not defined $the_seq_length){$the_seq_length = length ($result_hash{"N1A00"}{"seq"})};
		
		### BBBB
		### $the_seq_length
		
	};

	
};

### $the_seq_length


my %TOtal_mutation;
my %ID_times;

my $filespecies;
### $filespecies
### $tree_file
if ($tree_file =~ m/\w\w\-iT\_(.+)\_segment/xms){

	$filespecies = lc $1;

}elsif ($tree_file =~ m/\w\w\-iT\_*(.+)\_ALIGNED/xms){

	$filespecies = lc $1;
	### ok
}else{
	
	$filespecies = lc $tree_file;
	$filespecies =~ s/\.trees//g;
};

### $filespecies

if  ( not $Info_hash{$filespecies}){
	
	my @filespecies = split "_", $filespecies;
	## @filespecies
	my @tmp_infos = keys %Info_hash;
	## @tmp_infos
	for my $tmp_filespecies (@filespecies){
		

		my @tmp_grepage = grep (/${tmp_filespecies}/, @tmp_infos);

		@tmp_infos = @tmp_grepage;
		my $tmp_filespecies = "na";
	
		if ( scalar @tmp_infos  > 0 ){
			$tmp_filespecies =shift @tmp_infos;
		};
	
		if(defined($Info_hash{$tmp_filespecies})){
		
			$filespecies = $tmp_filespecies;
		
		}else{
			if ($filespecies){
				
				$Info_hash{$filespecies}{"Strain"}= "Ambigous";
				
			};
			
		};
	
	};
	
		
	
};


my $N_total = 0;
my $full_time_values=0;

# %result_hash


TIME:
for my $id (keys %result_hash){
	
	## $id 
	## X:$result_hash{$id}{"Date"} 
	$full_time_values = $full_time_values + $result_hash{$id}{"Date"} if $id ne "N0A0" ;
	my $racourcisseur = $id;
	## $id
	$racourcisseur =~ s/\N\d+A//;
	## $racourcisseur
	my $true_time = 0;
	
	POPOPO:	
	while ($racourcisseur){
					
		my $a = (length $racourcisseur) - 1;
		my $search_ID = join "", "N",$a,"A",$racourcisseur;
		
		## $search_ID
		


		
		my $time_value = $result_hash{$search_ID}{"Date"};
		
		## X : %{result_hash{$search_ID}}
		## $time_value
		my $true_time1 = $true_time;
		
		## $true_time1
		my $true_time2 = $true_time1 + $time_value;
		
		## $true_time2
		$true_time = $true_time2;		
		
		## $true_time
		$racourcisseur =~ s/.{1}$//;
		
				
		if ($search_ID =~ m/N1A/xms){ 
			
			$ID_times{$id}= $true_time;
			next TIME;
		
		};
			
	};
		
};


my $time_values = max values %ID_times;


# %result_hash

my $kj=0;

my @Measuring_timescale;
my @Measuring_timescale2;
for my $id (keys %ID_feuilles){
	$kj++;
	## $kj
	my $tmp_time=$ID_times{$id};
	push @Measuring_timescale,$tmp_time;
	my $tmp_time2 = $time_values - $ID_times{$id};
	push @Measuring_timescale2,$tmp_time2;
	
	## $id
	## X:$ID_feuilles{$id}
	## $tmp_time
	## $tmp_time2
	
};



my $Measuring_timescale = max(@Measuring_timescale) - min(@Measuring_timescale);
my $Measuring_timescale2 = max(@Measuring_timescale2) - min(@Measuring_timescale2);
##  @Measuring_timescale
##  @Measuring_timescale2
##  $Measuring_timescale
##  $Measuring_timescale2

my $species = "NA";
my $group= "NA";
my $genus= "NA";
my $family = "NA";	
				### $Info_hash{$filespecies}
if($Info_hash{$filespecies}){
					## AA				
	$group = $Info_hash{$filespecies}{"Group"};
	$genus= $Info_hash{$filespecies}{"Genus"};
	$family= $Info_hash{$filespecies}{"Family"};
	$species = $Info_hash{$filespecies}{"Strain"};					
					
};



my $seq_leni = length ($result_hash{"N1A00"}{"seq"}) + 1;
my %matrix_of_position;
tie %matrix_of_position, 'Tie::IxHash';

for (my $u=0; $u < $seq_leni;$u++){

	$matrix_of_position{$u}=0;

};

my $outfile_report = join "" , $outname,'_mut-report_CONSENSUS.txt';
open my $out, '>',  $outfile_report;
say {$out} join  "\t", "#ID","Position","Ref", "Mut", "Contexte","Type","Type2","Nature","Time","Time_full_branch","Time_tree","Seq_length", "Group","Genus","Family","Species","Sub";


my $outfile_fasta = join "" , $outname,'_seq.fasta';
open my $out4, '>',  $outfile_fasta;



#say {$out4} join  "", ">N1A00","\n",$result_hash{"N1A00"}{"seq"} ;
my $selected_count=0;

TAST:
for my $id (keys %result_hash){
	
	if( $id eq "N0A0" ){next TAST};
	my %res;
	my $time_data = 0;
	
	if( $id ne "N0A0" ){
	
		$time_data = $ID_times{$id};
	
	};
	
	## $id
	my $seq = $result_hash{$id}{"seq"};
	my $nat = $result_hash{$id}{"Nature"};
	say {$out4} join  "", ">",(join "|", $id,$group,$genus,$family,$species,$filespecies,$full_time_values,$time_data,$the_seq_length,$nat),"\n",$seq ;
	
	
	if ($result_hash{$id}{"anc_seq"}){
	my $trueseq=$seq;
	$trueseq=~ s/\-//g;
	my $realseqlength = $the_seq_length; 

	
	
	my $len = length $seq;

	
	## $tr_len_seq
	## $the_seq_length 
	
	$seq=~ s/N/-/g;
	## seq
	
	PRE:
	for(my $i = 0 ; $i <= $len ; $i += 1){
		
		my $last = 0;
		my $codon = substr ($seq, $i, 3);
		
		if ($codon =~  m/\-(\w\w)/xms){
			
			$codon=join "", $last, $1;
		
		};

		if ($codon =~  m/^\w/xms){
			
			if ($codon =~  m/\w(\w)\-/xms){
			
				$last = $1;
				next PRE; 
			
			}
						
			if ($codon =~  m/\w(\w)\w/xms){
					$res{$i}{$id} = {
						"Nucl" => $1,
						"Codon" => $codon,
					};
			}
				## $codon
				## $nuc
		}
		## $codon
		next PRE;
	}

## %res;

my $i = -1;

my @DNA_ref;
my $seq_rf  = $result_hash{$id}{"anc_seq"};

## TATAP:$id

@DNA_ref = split //, $seq_rf;

my $statut = $result_hash{$id}{"Nature"};
my $time_data = $time_values - $ID_times{$id};

## $seq_rf
## @DNA_ref

my $total_mut=0;



TURN:
while(@DNA_ref){
	my $nec = shift @DNA_ref;
	### $nec
	## $i
	
	if ($nec =~ m/\w/xms){ 
	for my $id (keys %{$res{$i}}){
				###$id
		
		if($res{$i}{$id}{"Nucl"} ne $nec){
					## x:$res{$i}{$id}{"Nucl"}
			my $up;
			my $down;
			for($res{$i}{$id}{"Codon"} =~ m/(\w)\w(\w)/xms){$up = $1; $down = $2};  
				
				
				my $u1 = $i;
				my $up1 = substr ($seq_rf, $u1, 1);
				UP1:
				while ($up1 =~ m/\-/xms ){
					$u1--;
					$up1 = substr ($seq_rf, $u1, 1);
				};
				
				my $u2 = $u1 - 1;
				my $up2 = substr ($seq_rf, $u2, 1);
				UP2:
				while ($up2 =~ m/\-/xms ){
					$u2--;
					$up2 = substr ($seq_rf, $u2, 1);
				};
				
				
				
				my $d1 = $i+2;
				my $down1 = substr ($seq_rf, $d1, 1);
				
				DOWN1:
				while ($down1 =~ m/\-/xms ){
					$d1++;
					$down1 = substr ($seq_rf, $d1, 1);
				};
				
				my $d2 = $d1+1;
				my $down2 = substr ($seq_rf, $d2, 1);
				
				DOWN2:
				while ($down2 =~ m/\-/xms ){
					$d2++;
					$down2 = substr ($seq_rf, $d2, 1);
				};
						 				

				
				## $group
				
				
				$N_total++;
				## $id
				
				say {$out} join "\t",  $id,$nec, $i+2, $res{$i}{$id}{"Nucl"}, $res{$i}{$id}{"Codon"}, (join "",$up2,$up1,"[",$nec,">",$res{$i}{$id}{"Nucl"},"]",$down1,$down2),
				(join "",$nec,">",$res{$i}{$id}{"Nucl"}),$statut,$time_data,$full_time_values,$time_values,$realseqlength, $group,$genus,$family,$species,$filespecies;
				
				my $corrected_i=$i+1;
				
				$matrix_of_position{$corrected_i}++;
				
				
				if ($id !~ m/N1A/xms){
				
					if($statut eq "Ancestor" ){
					
						$selected_count++;
					
					};
					
				};
				
				$TOtal_mutation{$id}{$corrected_i}={
					
					"Nuc-2"=>$up2,
					"Nuc-1"=>$up1,
					"up"=>$up,
					
					"Nuc+1"=>$down1,
					"Nuc+2"=>$down2,
					"down"=>$down,
					
					"Nuc"=>$nec,
					"Mut"=>$res{$i}{$id}{"Nucl"},
					"time"=>$time_data,
					"statut"=>$statut,
					"time"=>$time_data,
					"length"=>$realseqlength,
					"group"=>$group,
					"genus"=>$genus,
					"family"=>$family,
					"species"=>$species,
					"filespecies"=>$filespecies		
				
				
					
				};

			next;
				
			};
			
		};
			
	};
	
	$i = $i + 1;

	next TURN;
}

};
};

close $out;
close $out4;

#%TOtal_mutation





# ################################# Context Single ###################################################

open my $inTE, '<', $outfile_report;
my %matrix_of_singl;
tie %matrix_of_singl, 'Tie::IxHash';
TARN:		
while (my $line = <$inTE> ){
	chomp $line;
	my $mutation;
	
	if ($line =~ m/\w\[(\w\>\w)\]\w/xms){
		
		$mutation = $1;

		if ($matrix_of_singl{$mutation}) {$matrix_of_singl{$mutation}++}
		else {$matrix_of_singl{$mutation} = 1};
		
		next TARN;
		
	};
		
	next TARN;
}

close $inTE;



my $r = 'A,T,G,C';
my $kmersingle;
my @trikmersingle1;

for  my $kmer (glob"{$r}>{$r}"){
	
	for ($kmer =~ m/(\w)>(\w)/xms){$kmersingle ="$1>$2"};
	if ($1 ne $2){ push @trikmersingle1,$kmersingle };
	
}


my $outname_single = join "", $outname, "_matrix_single.txt";
open my $out_fin_single, '>',  $outname_single;

say {$out_fin_single} join "\t","Mutation Types", $outname;

my @trikmersingle2 = sort @trikmersingle1;



foreach my $tres(@trikmersingle2){
	## $tres
	if ($matrix_of_singl{$tres}) {
		## toto
		say {$out_fin_single} join "\t",$tres,$matrix_of_singl{$tres};
		}
	else {say {$out_fin_single} join "\t", $tres, 0};
};
	
close $out_fin_single;




my $outname_single_compil = join "", $outname, "_compilation_single.txt";
open my $out_fin_single_compil, '>',  $outname_single_compil;

say {$out_fin_single_compil} join "\t","ID","Mutation Types", "Value","Group","Genus","Family","Species","Full_age","Lenght","full_branch_time","tree_time";
foreach my $tres(@trikmersingle2){
	## $tres
	if ($matrix_of_singl{$tres}) {
		## toto
		say {$out_fin_single_compil} join "\t",$outname,$tres,$matrix_of_singl{$tres},$group,$genus,$family,$species,$the_seq_length,$full_time_values,$time_values;
		}
	else {say {$out_fin_single_compil} join "\t", $outname,$tres, 0,$group,$genus,$family,$species,$the_seq_length,$full_time_values,$time_values;
		};
};
	


close $out_fin_single_compil;

# ################################# Context Tri ###################################################


open my $inTEtri, '<', $outfile_report;
my %matrix;
tie %matrix, 'Tie::IxHash';

TARN:		
while (my $line = <$inTEtri> ){
		chomp $line;
	
	my $mutation;
	
	if ($line =~ m/(\w\[\w\>\w\]\w)/xms){
		
		$mutation = $1;
					

		if ($matrix{$mutation}) {$matrix{$mutation}++}
		else {$matrix{$mutation} = 1};
		
		next TARN;
	};
		
		next TARN;
}

close $inTEtri;


my $kmero;
my @trikmer1;
for  my $kmer (glob"{$r}{$r}>{$r}{$r}"){
	## $kmer
	for ($kmer =~ m/(\w)(\w)>(\w)(\w)/xms){$kmero =join "",$1,'[',$2,'>',$3,']',$4};
	if ($2 ne $3){ push @trikmer1,$kmero };
}


my $outname_final = join "", $outname, "_matrix_line.txt";
open my $out_fin, '>',  $outname_final;

say {$out_fin} join "\t","Mutation Types", $outname;

## @trikmer1;

my @trikmer2 = sort @trikmer1;



foreach my $tres(@trikmer2){
	## $tres
	if ($matrix{$tres}) {
		## toto
		say {$out_fin} join "\t",$tres,$matrix{$tres};
		}
	else {say {$out_fin} join "\t", $tres, 0};
};

close $out_fin;


# ##### Read matrix96


open my $inTE2, '<', $outfile_report;
my %matrix2;
tie %matrix2, 'Tie::IxHash';


TARN:		
while (my $line = <$inTE2> ){
		chomp $line;
	
	my $mutation;

	if ($line =~ m/(\w\[\w\>\w\]\w)/xms){
		
		$mutation = $1;

		if ($mutation =~ m/G>T|G>C|G>A|A>T|A>G|A>C/xms){
			
			if ($line =~ m/(\w)(\[\w\>\w\])(\w)/xms){
			
				$mutation = join "" ,$3,$2,$1;
			
			};
			
			$mutation =~ tr/ACGTacgt/TGCAtgca/;
		
		};
			
					

		if ($matrix2{$mutation}) {$matrix2{$mutation}++}
		else {$matrix2{$mutation} = 1};
		
		next TARN;
	};
		
		next TARN;
}



close $inTE2;


$r = 'A,T,G,C';
my $t = 'C';
my $kmeroNA;
my @trikmerNA;
for  my $kmer (glob"{$r}{$t}>{$r}{$r}"){
	for ($kmer =~ m/(\w)(\w)>(\w)(\w)/xms){$kmeroNA = join "",$1,'[',$2,'>',$3,']',$4};
	if ($2 ne $3){ push @trikmerNA,$kmeroNA };
}

@trikmer2 = sort @trikmerNA;

$t = 'T';
my $kmeroNA2;
my @trikmerNA3;
for  my $kmer (glob"{$r}{$t}>{$r}{$r}"){
	for ($kmer =~ m/(\w)(\w)>(\w)(\w)/xms){$kmeroNA2 =join "",$1,'[',$2,'>',$3,']',$4};
	if ($2 ne $3){ push @trikmerNA3,$kmeroNA2 };
}

my @trikmer4 = sort @trikmerNA3;

my @trikmer5 = (@trikmer2, @trikmer4);


my $outname_final3 = join "", $outname, "_matrix96_line.txt";
open my $out_fin3, '>',  $outname_final3;

say {$out_fin3} join "\t","Mutation Types", $outname;


foreach my $tres(@trikmer5){

	if ($matrix2{$tres}) {

		say {$out_fin3} join "\t",$tres,$matrix2{$tres};
		}
	else {say {$out_fin3} join "\t", $tres, 0};
};

close $out_fin3;

# ################################# Context Penta ###################################################




open my $inTE_penta, '<', $outfile_report;
my %matrix_penta;
tie %matrix_penta, 'Tie::IxHash';
my $kmeroPenta ;
for  my $kmer (glob"{$r}{$r}{$r}>{$r}{$r}{$r}"){
	
	for ($kmer =~ m/(\w)(\w)(\w)\>(\w)(\w)(\w)/xms){$kmeroPenta =join "", $1,$2,'[',$3,'>',$4,']',$5,$6};
	if ($3 ne $4){ $matrix_penta{$kmeroPenta} = 0 };
}


TARN:		
while (my $line = <$inTE_penta> ){
		chomp $line;
	
	my $mutation;
	
	if ($line =~ m/(\w\w\[\w\>\w\]\w\w)/xms){
		
		$mutation = $1;

		$matrix_penta{$mutation}++;
		
		next TARN;
	};
		
		next TARN;
}

close $inTE_penta;

## %matrix_penta


my $outname_final_penta = join "", $outname, "_matrix_penta_line.txt";
open my $out_fin_penta, '>',  $outname_final_penta;

say {$out_fin_penta} join "\t","Mutation Types", $outname;


for my $tres(keys %matrix_penta){

		say {$out_fin_penta} join "\t",$tres,$matrix_penta{$tres};

};


close $out_fin_penta;





# ################################    Search doubles   ##############################################



open my $inTE3, '<', $outfile_report;

my %di_mut_matrix;
tie %di_mut_matrix, 'Tie::IxHash';

my $outfile_report_double = join "" , $outname,'_mut-report_DOUBLE.txt';
open my $out_double, '>',  $outfile_report_double;
say {$out_double} join  "\t", "#ID","Position","Ref", "Mut","Type","Type2","Nature","Time","Seq_length", "Group","Genus","Family","Species","Sub";


for my $od (keys %TOtal_mutation){

	for my $posi (keys %{$TOtal_mutation{$od}}){
	
		my $last = $posi-1;
		my $next = $posi+1;
		
		if ($TOtal_mutation{$od}{$last}){
				
				
				my $L1=$TOtal_mutation{$od}{$last}{"Nuc-2"};
				my $L2=$TOtal_mutation{$od}{$last}{"Nuc-1"};
				my $R1=$TOtal_mutation{$od}{$last}{"Nuc"};
				my $R2=$TOtal_mutation{$od}{$posi}{"Nuc"};
				my $M1=$TOtal_mutation{$od}{$last}{"Mut"};
				my $M2=$TOtal_mutation{$od}{$posi}{"Mut"};
				my $N1=$TOtal_mutation{$od}{$posi}{"Nuc+1"};
				my $N2=$TOtal_mutation{$od}{$posi}{"Nuc+2"};
				my $statut_double =$TOtal_mutation{$od}{$posi}{"statut"};
				my $time_double =$TOtal_mutation{$od}{$posi}{"time"};
				my $realseqlength_double =$TOtal_mutation{$od}{$posi}{"length"};
				my $group_double =$TOtal_mutation{$od}{$posi}{"group"};
				my $genus_double =$TOtal_mutation{$od}{$posi}{"genus"};
				my $family_double =$TOtal_mutation{$od}{$posi}{"family"};
				my $species_double =$TOtal_mutation{$od}{$posi}{"species"};
				my $filespecies_double = $TOtal_mutation{$od}{$posi}{"filespecies"};
				
				
				my $Dim=join "",$R1,$R2,">",$M1,$M2; 
				my $DiMut=join "",$L2,"[",$Dim,"]","$N1"; 
				my $DiMut_penta=join "", $L1,$L2,"[",$R1,$R2,">",$M1,$M2,"]","$N1","$N2"; 
				
				if ($di_mut_matrix{$Dim}) {$di_mut_matrix{$Dim}++}
				else {$di_mut_matrix{$Dim} = 1};
				
				say {$out_double} join "\t",  $od,$posi, (join "" ,$R1,$R2),(join "",$M1,$M2), $DiMut_penta,$Dim,
				$statut_double,$time_double,$realseqlength_double, $group_double,$genus_double,$family_double,$species_double,$filespecies_double;

	
		}
	
	}

}

close $out_double;


my $kmero_double;
my @trikmer_double;

for  my $kmer (glob"{$r}{$r}>{$r}{$r}"){

	my $dimer_out;
	my $dimer_in;
	
	for ($kmer =~ m/(\w)(\w)>(\w)(\w)/xms){
		$dimer_in=join "" ,$1,$2; 
		$dimer_out=join "" ,$3,$4; 
		$kmero_double ="$dimer_in>$dimer_out";
	};
	
	if ($dimer_in ne $dimer_out){ push @trikmer_double,$kmero_double };

}

@trikmer2 = sort @trikmer_double;



my $outname_final_double3 = join "", $outname, "_matrix_double_line.txt";
open my $out_fin_double2, '>',  $outname_final_double3;

say {$out_fin_double2} join "\t","Mutation Types", $outname;



foreach my $tres(@trikmer2){
	## $tres
	if ($di_mut_matrix{$tres}) {
		## toto
		say {$out_fin_double2} join "\t",$tres,$di_mut_matrix{$tres};
		}
	else {say {$out_fin_double2} join "\t", $tres, 0};
};

close $out_fin_double2;

# ################################    Search PingPong   ##############################################

my %partition;
tie %partition, 'Tie::IxHash';

my $kmero_full;
my %trikrmer_full_PING_PONG;
tie %trikrmer_full_PING_PONG, 'Tie::IxHash';

my %trikrmer_full_NOPING_PONG;
tie %trikrmer_full_NOPING_PONG, 'Tie::IxHash';

my %Ping_sum;
my %Ping_sum_no_rest;



my %Ping_pos_sum;


LECTURE:
for my $id (keys %ID_feuilles){
	
	my $racourcisseur = $id;
	$racourcisseur =~ s/N\d+A//;
	
	my $seqId = $ID_feuilles{$id};
	
	while ($racourcisseur){
			
			my $a= (length $racourcisseur) - 1;
			my $search_ID = join "", "N",$a,"A",$racourcisseur;
			my $mark = 0;

			POSITION:
			for my $pos (keys %{$TOtal_mutation{$search_ID}}){
				
				
				## $search_ID
				
				$partition{$seqId}{$pos}{$search_ID}=$TOtal_mutation{$search_ID}{$pos};
					 
				$mark = 1;
				
			};
		
		$racourcisseur =~ s/.{1}$//;
				
	};

};

my %Ping_infox;
	
for my $IDres (keys %partition){

	for my $position (keys %{$partition{$IDres}}){

		my @values = keys %{$partition{$IDres}{$position}};

		if (scalar @values > 1){
			
			@values= sort  @values;
			my $testpo=0;
			for my $value (@values){
			
				my $tri_kmero_id=join "",$partition{$IDres}{$position}{$value}{"Nuc-1"},"[",$partition{$IDres}{$position}{$value}{"Nuc"},
				">",$partition{$IDres}{$position}{$value}{"Mut"},"]",$partition{$IDres}{$position}{$value}{"Nuc+1"};
				
				my $statut = $result_hash{$value}{"Nature"};
				
				my $ID_TmP = join "_", $value,$position,$statut;
				my $ID_TmP2 = join "_", $value,$position,$IDres;
				
				if ($partition{$IDres}{$position}{$values[1]}{"Mut"} eq $partition{$IDres}{$position}{$values[0]}{"Nuc"}){
				
					
					$Ping_sum{$ID_TmP}=1 ;

					my $tri_kmero_id=join "",$partition{$IDres}{$position}{$value}{"Nuc-1"},"[",$partition{$IDres}{$position}{$value}{"Nuc"},
					">",$partition{$IDres}{$position}{$value}{"Mut"},"]",$partition{$IDres}{$position}{$value}{"Nuc+1"};
					## POP
					$trikrmer_full_PING_PONG{$tri_kmero_id}{$ID_TmP}=0;
					$Ping_pos_sum{$position}{$ID_TmP}=0;

				}else{
						
					$Ping_sum_no_rest{$ID_TmP}=1;
					$trikrmer_full_NOPING_PONG{$tri_kmero_id}{$ID_TmP}=0;
				
				};
				
				my $testpopi =$testpo +1;
				
				if ($values[$testpopi]){
					
				
					my $ID1 = $values[$testpo];
					my $ID2 = $values[$testpopi];
				
					my $time1 = $TOtal_mutation{$ID1}{$position}{"time"};
					my $time2 = $TOtal_mutation{$ID2}{$position}{"time"};
					
					my $PING_score = $time1 +  $time2;
					
					## $PING_score
					
					$Ping_infox{$ID_TmP2}={
						"Branche" => $IDres,
						"Time" =>  $PING_score,
						"Position" => $position
						
					};
					
				};
				
				$testpo++;
				
			};	
			
							
		}else{
			
			my $value = shift @values;
			my $ID_TmP = join "_", $value,$position;
			my $tri_kmero_id=join "",$partition{$IDres}{$position}{$value}{"Nuc-1"},"[",$partition{$IDres}{$position}{$value}{"Nuc"},
			">",$partition{$IDres}{$position}{$value}{"Mut"},"]",$partition{$IDres}{$position}{$value}{"Nuc+1"};
				
			$trikrmer_full_NOPING_PONG{$tri_kmero_id}{$ID_TmP}=0;
				
		};
		
	};
	
};


## %trikrmer_full_NOPING_PONG

my $Ping_sum = 0;
if (%Ping_sum){
	$Ping_sum = sum values %Ping_sum;
};


my $Ping_sum_select = 0;

if (%Ping_sum){
	
	for my $keyss ( keys %Ping_sum ){

		if ($keyss !~ m/N1A/xms){
			
			if ($keyss =~ m/Ancestor/xms){
					
				$Ping_sum_select ++;
					
			};
					
		};
	
	};

};



my $Ping_sum_no_rest = 0;
if (%Ping_sum_no_rest){
	$Ping_sum_no_rest = sum values %Ping_sum_no_rest;
};

my $Ratio_PING=0;

if($Ping_sum > 0){
	
	$Ratio_PING = (sprintf  "%.3f", ($Ping_sum / $N_total))*100 ;	

};



my $r3 = 'A,T,G,C';
my $kmero3 ;
my @trikmersingle3;

for  my $kmer (glob"{$r3}{$r3}>{$r3}{$r3}"){
	
	if ($kmer =~ m/(\w)(\w)\>(\w)(\w)/xms){
			
		$kmero3 =join "", $1,'[',$2,'>',$3,']',$4;
		
		if ($2 ne $3){ 
			
			push @trikmersingle3,$kmero3;
			
		};
	
	};

};  




###########

#my $recombo_file = join "", $outname, "_Recombo_pose.txt";
#open my $out_recombo, '>',  $recombo_file;
#say {$out_recombo} join "\t","PAIRE_ID","DIST_TIME","DIST_POSE,Branche";


#foreach my $Paire(keys %Ping_infox){
	
#	my $local_pose = $Ping_infox{$Paire}{"Position"};
#	my $local_time = $Ping_infox{$Paire}{"Time"};
#	my $locale_branch = $Ping_infox{$Paire}{"Branche"};
	
#	foreach my $Paire2(keys %Ping_infox){
		
#		if ($Ping_infox{$Paire2}{"Branche"}==$locale_branch){
			
#			my $local_pose2 = $Ping_infox{$Paire2}{"Position"};
#			my $local_time2 = $Ping_infox{$Paire2}{"Time"};
			
#			my $Dist_time = abs($local_time -$local_time2);	
#			my $Dist_pos = abs($local_pose -$local_pose2);	
			
#			my $comp_id = join "_", $Paire , $Paire2;
			
#			say {$out_recombo} join "\t",$comp_id,  $Dist_time, $Dist_pos, $locale_branch ;
		
#		};
		
#	};
	
#};

#############

###########

#my $recombo_file2 = join "", $outname, "_Recombo_pose2.txt";
#open my $out_recombo2, '>',  $recombo_file2;
#say {$out_recombo2} join "\t","PAIRE_ID","DIST_TIME","DIST_POSE,Branche";


#foreach my $Paire(keys %Ping_infox){
	
#	my $local_pose = $Ping_infox{$Paire}{"Position"};
#	my $local_time = $Ping_infox{$Paire}{"Time"};
#	my $locale_branch = $Ping_infox{$Paire}{"Branche"};
	
#	say {$out_recombo2} join "\t",$Paire,  $local_time, $local_pose, $locale_branch ;
		

#};

#############




my @trikmersingle4 = sort @trikmersingle3;


my $outname_tri192_NOping_full = join "", $outname, "_matrix_192_NO_PING-PONG_full.txt";
open my $out_fin_tri_NOping192_full, '>',  $outname_tri192_NOping_full;
say {$out_fin_tri_NOping192_full} join "\t","Mutation Types", $outname;


foreach my $NOPing_keys_192_full(@trikmersingle4){
		
		if ($trikrmer_full_NOPING_PONG{$NOPing_keys_192_full}) {

			say {$out_fin_tri_NOping192_full} join "\t",$NOPing_keys_192_full,  scalar values %{$trikrmer_full_NOPING_PONG{$NOPing_keys_192_full}} ;		
			
		}else {
			
			say {$out_fin_tri_NOping192_full} join "\t",$NOPing_keys_192_full,  0 ;

		};
		
};


close $out_fin_tri_NOping192_full;




my $outname_tri192_ping_full = join "", $outname, "_matrix_192_PING-PONG_full.txt";
open my $out_fin_tri_ping192_full, '>',  $outname_tri192_ping_full;
say {$out_fin_tri_ping192_full} join "\t","Mutation Types", $outname;


foreach my $Ping_keys_192_full(@trikmersingle4){
		
		if ($trikrmer_full_PING_PONG{$Ping_keys_192_full}) {

			say {$out_fin_tri_ping192_full} join "\t",$Ping_keys_192_full,  scalar values %{$trikrmer_full_PING_PONG{$Ping_keys_192_full}} ;		
			
		}else {
			
			say {$out_fin_tri_ping192_full} join "\t",$Ping_keys_192_full,  0 ;

		};
		
};

close $out_fin_tri_ping192_full;











#############################################  position rate  #####################


my $outname_posrate = join "", $outname, "_position_rate.txt";
open my $out_posrate, '>',  $outname_posrate;

#say {$out_posrate} join "\t","ID","Position","nMutation", "Rate","full_branch_time","tree_time","Group","Genus","Family","Species";
say {$out_posrate} join "\t","#Position",$outname;

for (my $u=1; $u < $seq_leni;$u++){

	my $popol = $matrix_of_position{$u};
	my $popol_rate = $popol/ $full_time_values;
	#say {$out_posrate} join "\t",$outname,$u,$popol,$popol_rate,$full_time_values,$time_values,$group,$genus,$family,$species;
	say {$out_posrate} join "\t",$u,$popol_rate;
};

close $out_posrate;



my $outname_POS_ping_full = join "", $outname, "_matrix_POSITION_PING-PONG_full.txt";
open my $out_fin_POS_full, '>',  $outname_POS_ping_full;
say {$out_fin_POS_full} join "\t","#id", 'Position',"TT","PING";

## %Ping_pos_sum
TOUR:
for (my $u=1; $u < $seq_leni;$u++){
	
	if ( $matrix_of_position{$u} == 0 ){ next TOUR };
				
	print {$out_fin_POS_full} join "\t",$outname, $u,   $matrix_of_position{$u} ;		
		
	print {$out_fin_POS_full} "\t";
		
	if (%Ping_pos_sum) {

		if ($Ping_pos_sum{$u}) {
				
			say {$out_fin_POS_full} scalar values %{$Ping_pos_sum{$u}} ;		
		
		}else{
			
			say {$out_fin_POS_full} "0" ;

		};
	
	}else {
			
		say {$out_fin_POS_full} "0" ;

	};

};

close $out_fin_POS_full;
%matrix_of_position=();
# ################## Analyse Sub Rev Evo rates mono context ###############################


# ################ matrix compo 


open my $in_compo_mono, '<', $outfile_fasta;
my %seq_for_compo_mono;
tie %seq_for_compo_mono, 'Tie::IxHash';
my $seq_id_compo_mono;
my $seq_compo_mono;

LINE:		
while (my $line = <$in_compo_mono> ){	
	chomp $line;
	if ($line =~ m/>(.*)/xms) {
			
		 # add current seq to hash  (if any)
		if($seq_compo_mono) {
			$seq_for_compo_mono{$seq_id_compo_mono} = $seq_compo_mono;	
			$seq_compo_mono = q{};
		}
		
		# extract new seq_id
		$seq_id_compo_mono = $1;
		next LINE;
			
	}
	# elongate current seq (seq can be broken sevetal lines)
	$seq_compo_mono .= $line;
}

#add last seq to hash (if any)
$seq_for_compo_mono{$seq_id_compo_mono} = $seq_compo_mono if $seq_compo_mono;



close $in_compo_mono;





	
my $k = 'A,T,G,C';
my $kmero_compo_mono;

my %composition_tri_mono;
tie %composition_tri_mono, 'Tie::IxHash';



for  my $kmer_compo_mono (glob"{$k}"){
	$composition_tri_mono{$kmer_compo_mono}=0;
}

my $full_mono_compo  = 0;


for my $id_compo_mono (keys %seq_for_compo_mono){

	$seq_compo_mono = $seq_for_compo_mono{$id_compo_mono};
	
	for my $kmer_compo_mono (keys %composition_tri_mono){
		
		my $countage_compo_mono = () = $seq_compo_mono =~ /$kmer_compo_mono/g;

		$composition_tri_mono{$kmer_compo_mono}+= $countage_compo_mono;
		$full_mono_compo +=  $countage_compo_mono;
		
	}
}



my $out_matrix_compo_mono = join "" , $outname,'_matrix_compo_mono.txt';
open my $out_compo_mono, '>',  $out_matrix_compo_mono;
say {$out_compo_mono} join  "\t", "#Trik","Count";


for my $last_compting_mono (keys %composition_tri_mono){
	
	say {$out_compo_mono} join  "\t", $last_compting_mono,$composition_tri_mono{$last_compting_mono};
	
};

close $out_compo_mono;



my %mix_matrix_mono;
tie %mix_matrix_mono, 'Tie::IxHash';			
open my $in_rec_mono, '<', $outname_single;


my $line_mono = <$in_rec_mono> ;


LINE_mono:
while (my $line_mono = <$in_rec_mono> ){

	chomp $line_mono;
	if ($line_mono =~ m/Mutation\ Types/xms){next LINE_mono};
	
	## $line
	
	my @values = split "\t" ,$line_mono;
	my $remove = shift @values;
	my $mediane= shift @values;
	
	## $remove
	## $mediane
	
	$mix_matrix_mono{$remove} = $mediane;
		
};
		
close $in_rec_mono;



## %mix_matrix

my $out_matrix_compo_tx_mono = join "" , $outname,'_matrix_txS_mono.txt';
open my $out_compo_tx_mono, '>',  $out_matrix_compo_tx_mono;
say {$out_compo_tx_mono} join  "\t", "Mutation Types",$outname;


my %TxS_mono;
tie %TxS_mono, 'Tie::IxHash';
my $TxS_mono=0;


for my $id_rates_mono (keys %mix_matrix_mono){
	
	my $count_value = $mix_matrix_mono{$id_rates_mono};
	my $tri_compo_kmer;
	if ($id_rates_mono =~ m/(\w)\>\w/xms){$tri_compo_kmer=join "", $1};
	my $value_of_compo = $composition_tri_mono{$tri_compo_kmer};
	
	#my $mut_tx = $count_value / ( $value_of_compo * $full_time_values);
	my $mut_tx;
	my $diviseur = $value_of_compo * $full_time_values;
	
	if ($count_value > 0 && $diviseur > 0 ){
		$mut_tx = $count_value / $diviseur;
	}else{
		$mut_tx = 0;
	};


	say {$out_compo_tx_mono} join  "\t", $id_rates_mono , $mut_tx;
	$TxS_mono{$id_rates_mono}=$mut_tx;
	$TxS_mono = $TxS_mono + $mut_tx;
	
	my $prop_compo = $value_of_compo/$full_mono_compo;
	my $Txs_corr = $mut_tx *  $prop_compo;
	$TxS_mono= $TxS_mono + $Txs_corr;
	
};


close $out_compo_tx_mono;






my $outname_single_compil_txs = join "", $outname, "_compilation_TxS_mono.txt";
open my $out_fin_single_compil_txs, '>',  $outname_single_compil_txs;




say {$out_fin_single_compil_txs} join "\t","ID","Lifs","Mutation Types", "Value","Group","Genus","Family","Species","Full_age","Lenght","full_branch_time","tree_time";

for my $id_TxS_mono (keys %TxS_mono){
	
	my $val_TxS = $TxS_mono{$id_TxS_mono};
	
	say {$out_fin_single_compil_txs} join "\t",$outname,$id_TxS_mono,$val_TxS,$group,$genus,$family,$species,$the_seq_length,$full_time_values,$time_values;

};
	


close $out_fin_single_compil_txs;












my $out_matrix_compo_txR_mono = join "" , $outname,'_matrix_txR_mono.txt';
open my $out_compo_txR_mono, '>',  $out_matrix_compo_txR_mono;
say {$out_compo_txR_mono} join  "\t", "Mutation Types",$outname;


my $out_matrix_compo_txE_mono = join "" , $outname,'_matrix_txE_mono.txt';
open my $out_compo_txE_mono, '>',  $out_matrix_compo_txE_mono;
say {$out_compo_txE_mono} join  "\t", "Mutation Types",$outname;



my $group_rate = $Info_hash{$filespecies}{"Group"};
my $genus_rate = $Info_hash{$filespecies}{"Genus"};
my $family_rate = $Info_hash{$filespecies}{"Family"};
my $species_rate = $Info_hash{$filespecies}{"Strain"};	
					
					


my %TxR_mono;
tie %TxR_mono, 'Tie::IxHash';
my %TxE_mono;
tie %TxE_mono, 'Tie::IxHash';
my $TxR_mono=0;
my $TxE_mono=0;

for my $id_TxS_mono (keys %TxS_mono){
	my $tri_compo_kmer;
	if ($id_TxS_mono =~ m/(\w)\>\w/xms){$tri_compo_kmer=join "", $1};
	my $value_of_compo = $composition_tri_mono{$tri_compo_kmer};
	
	my $val_TxS = $TxS_mono{$id_TxS_mono};
	my $id_TxR_mono;
	if ($id_TxS_mono =~ m/(\w)\>(\w)/xms){$id_TxR_mono=join "",$2,">",$1};
	my $val_TxR = $TxS_mono{$id_TxR_mono};
	
	my $true_revert;
	if ($val_TxS > $val_TxR){$true_revert=$val_TxR}else{$true_revert=$val_TxS};
	$TxR_mono{$id_TxS_mono}=$true_revert;
	my $val_TxE=$val_TxS-$true_revert;
	$TxE_mono{$id_TxS_mono}=$val_TxE;
	say {$out_compo_txR_mono} join  "\t", $id_TxS_mono , $true_revert;
	say {$out_compo_txE_mono} join  "\t", $id_TxS_mono , $val_TxE;
	
	
	
	my $prop_compo = $value_of_compo/$full_mono_compo;
	my $Txr_corr = $true_revert *  $prop_compo;
	$TxR_mono= $TxR_mono + $Txr_corr;
	my $Txe_corr = $val_TxE *  $prop_compo;
	$TxE_mono= $TxE_mono + $Txe_corr;
	
	
}


close $out_compo_txE_mono;
close $out_compo_txR_mono;
























# ################## Analyse Sub Rev Evo rates tri context ###############################


# ################ matrix compo 


open my $in_compo, '<', $outfile_fasta;
my %seq_for_compo;
tie %seq_for_compo, 'Tie::IxHash';
my $seq_id_compo;
my $seq_compo;
my $N_seq=0;

LINE:		
while (my $line = <$in_compo> ){	
	chomp $line;
	if ($line =~ m/>(.*)/xms) {
			
		 # add current seq to hash  (if any)
		if($seq_compo) {
			$seq_for_compo{$seq_id_compo} = $seq_compo;	
			$seq_compo = q{};
		}
		
		# extract new seq_id
		$seq_id_compo = $1;
		$N_seq++;
		next LINE;
			
	}
	# elongate current seq (seq can be broken sevetal lines)
	$seq_compo .= $line;
}

#add last seq to hash (if any)
$seq_for_compo{$seq_id_compo} = $seq_compo if $seq_compo;



close $in_compo;





	
$k = 'A,T,G,C';
my $kmero_compo;

my %composition_tri;
tie %composition_tri, 'Tie::IxHash';


for  my $kmer_compo (glob"{$k}{$k}{$k}"){
	#$composition_tri{$kmer_compo}=0;
	push @{$composition_tri{$kmer_compo}},0;
}


my @full_tri_compo;
my $number_tips = 1;
#my $number_tips = scalar keys %seq_for_compo;

for my $id_compo (keys %seq_for_compo){

	$seq_compo = $seq_for_compo{$id_compo};
	my $tmp_full_tri_compo =0;
	
	for my $kmer_compo (keys %composition_tri){
		
		my $countage_compo = () = $seq_compo =~ /$kmer_compo/g;
		
		push @{$composition_tri{$kmer_compo}},$countage_compo;
		$tmp_full_tri_compo +=  $countage_compo;
		
	}
	
	push @full_tri_compo,$tmp_full_tri_compo;
	
}

my $full_tri_compo = mean(@full_tri_compo);

my $out_matrix_compo = join "" , $outname,'_matrix_compo.txt';
open my $out_compo, '>',  $out_matrix_compo;
say {$out_compo} join  "\t", "#Trik","Count";


for my $last_compting (keys %composition_tri){
	
	say {$out_compo} join  "\t", $last_compting,mean(@{$composition_tri{$last_compting}});
	
};

close $out_compo;



my %mix_matrix;
tie %mix_matrix, 'Tie::IxHash';			
open my $in_rec, '<', $outname_final;


my $line = <$in_rec> ;


LINE:
while (my $line = <$in_rec> ){

	chomp $line;
	if ($line =~ m/Mutation\ Types/xms){next LINE};
	
	## $line
	
	my @values = split "\t" ,$line;
	my $remove = shift @values;
	my $mediane= shift @values;
	
	## $remove
	## $mediane
	
	$mix_matrix{$remove} = $mediane;
		
};
		
close $in_rec;



## %mix_matrix

my $out_matrix_compo_tx = join "" , $outname,'_matrix_txS.txt';
open my $out_compo_tx, '>',  $out_matrix_compo_tx;
say {$out_compo_tx} join  "\t", "Mutation Types",$outname;


my %TxS;
tie %TxS, 'Tie::IxHash';
my $total_TXS=0;

for my $id_rates (keys %mix_matrix){
	
	my $count_value = $mix_matrix{$id_rates};
	my $tri_compo_kmer;
	if ($id_rates =~ m/(\w)\[(\w)\>\w\](\w)/xms){$tri_compo_kmer=join "", $1,$2,$3};
	my $value_of_compo = mean(@{$composition_tri{$tri_compo_kmer}});
	
	#my $mut_tx = $count_value / ( $value_of_compo * $full_time_values);
	my $mut_tx;
	my $diviseur = ($value_of_compo / $number_tips) * $full_time_values;
	
	if ($count_value > 0 && $diviseur > 0 ){
		$mut_tx = $count_value / $diviseur;
	}else{
		$mut_tx = 0;
	};


	say {$out_compo_tx} join  "\t", $id_rates , $mut_tx;
	$TxS{$id_rates}=$mut_tx;
	
	$total_TXS += $count_value;
};


close $out_compo_tx;

	
	
	
	
	

my $out_matrix_compo_txR = join "" , $outname,'_matrix_txR.txt';
open my $out_compo_txR, '>',  $out_matrix_compo_txR;
say {$out_compo_txR} join  "\t", "Mutation Types",$outname;


my $out_matrix_compo_txE = join "" , $outname,'_matrix_txE.txt';
open my $out_compo_txE, '>',  $out_matrix_compo_txE;
say {$out_compo_txE} join  "\t", "Mutation Types",$outname;

	
					


my %TxR;
tie %TxR, 'Tie::IxHash';
my %TxE;
tie %TxE, 'Tie::IxHash';

my $total_TXR = 0;
my $total_TXE = 0;

for my $id_TxS (keys %mix_matrix){
	
	my $tri_compo_kmer;
	if ($id_TxS =~ m/(\w)\[(\w)\>\w\](\w)/xms){$tri_compo_kmer=join "", $1,$2,$3};
	my $value_of_compo = mean(@{$composition_tri{$tri_compo_kmer}});
	
	
	my $val_TxS = $mix_matrix{$id_TxS};
	my $id_TxR;
	if ($id_TxS =~ m/(\w)\[(\w)\>(\w)\](\w)/xms){$id_TxR=join "", $1,"[",$3,">",$2,"]",$4};
	my $val_TxR = $mix_matrix{$id_TxR};
	
	my $true_revert;
	my $Rate_true_revert=0;
	if ($val_TxS > $val_TxR){$true_revert=$val_TxR}else{$true_revert=$val_TxS};
	
	my $diviseur = ($value_of_compo / $number_tips) * $full_time_values;
	

	if ($val_TxR > 0 && $diviseur > 0 ){
		$Rate_true_revert = $true_revert / $diviseur;
	};

	$TxR{$id_TxS}=$Rate_true_revert;
	
	
	my $val_TxE=$val_TxS-$true_revert;	
	my $Rate_TxE=0;

	
	if ($val_TxE > 0 && $diviseur > 0 ){
		$Rate_TxE = $val_TxE / $diviseur;
	};

	$TxE{$id_TxS}=$Rate_TxE;

	say {$out_compo_txR} join  "\t", $id_TxS , $Rate_true_revert;
	say {$out_compo_txE} join  "\t", $id_TxS , $Rate_TxE;
	
	$total_TXR += $true_revert;
	$total_TXE += $val_TxE;


}

close $out_compo_txE;
close $out_compo_txR;


my $TxS_total = $total_TXS/ ($full_time_values * $the_seq_length);
my $TxR_total = $total_TXR/ ($full_time_values * $the_seq_length);
my $TxE_total = $total_TXE/ ($full_time_values * $the_seq_length);


# ################# Analyse Sub Rev Evo rates penta context ###############################
# ################ matrix compo 




open my $in_compo_penta, '<', $outfile_fasta;
my %seq_for_compo_penta;
tie %seq_for_compo_penta, 'Tie::IxHash';
my $seq_id_compo_penta;
my $seq_compo_penta;

LINE:		
while (my $line = <$in_compo_penta> ){	
	chomp $line;
	if ($line =~ m/>(.*)/xms) {
			
		 # add current seq to hash  (if any)
		if($seq_compo_penta) {
			$seq_for_compo_penta{$seq_id_compo_penta} = $seq_compo_penta;	
			$seq_compo_penta = q{};
		}
		
		# extract new seq_id
		$seq_id_compo_penta = $1;
		next LINE;
			
	}
	# elongate current seq (seq can be broken sevetal lines)
	$seq_compo_penta .= $line;
}

#add last seq to hash (if any)
$seq_for_compo_penta{$seq_id_compo_penta} = $seq_compo_penta if $seq_compo_penta;




close $in_compo_penta;





my $kmero_compo_penta;

my %composition_tri_penta;
tie %composition_tri_penta, 'Tie::IxHash';


for  my $kmer_compo (glob"{$k}{$k}{$k}{$k}{$k}"){
	$composition_tri_penta{$kmer_compo}=0;
}

my $full_penta_compo = 0;

for my $id_compo (keys %seq_for_compo_penta){

	$seq_compo_penta = $seq_for_compo_penta{$id_compo};
	
	for my $kmer_compo (keys %composition_tri_penta){
		
		my $countage_compo = () = $seq_compo_penta =~ /$kmer_compo/g;
		
		$composition_tri_penta{$kmer_compo}+= $countage_compo;
		$full_penta_compo +=  $countage_compo;
	}
}




my $out_matrix_compo_penta = join "" , $outname,'_matrix_compo_penta.txt';
open my $out_compo_penta, '>',  $out_matrix_compo_penta;
say {$out_compo_penta} join  "\t", "#Trik","Count";


for my $last_compting (keys %composition_tri_penta){
	
	say {$out_compo_penta} join  "\t", $last_compting,$composition_tri_penta{$last_compting};
	
};

close $out_compo_penta;





my %mix_matrix_penta;
tie %mix_matrix_penta, 'Tie::IxHash';			
open my $inrecpenta, '<', $outname_final_penta;


LINO:
while (my $lino = <$inrecpenta> ){
	
	chomp $lino;
	if ($lino =~ m/Mutation\ Types/xms){next LINO};
	

	my @values = split "\t" ,$lino;
	my $remove = shift @values;
	my $mediane= shift @values;

	
	$mix_matrix_penta{$remove} = $mediane;
	next LINO;
};
		
close $inrecpenta;












my $out_matrix_compo_penta_tx = join "" , $outname,'_matrix_txS_penta.txt';
open my $out_compo_penta_tx, '>',  $out_matrix_compo_penta_tx;
say {$out_compo_penta_tx} join  "\t", "Mutation Types",$outname;


my %TxS_penta;
tie %TxS_penta, 'Tie::IxHash';
my $total_TXS_penta=0;

for my $id_rates (keys %mix_matrix_penta){
	
	my $count_value = $mix_matrix_penta{$id_rates};
	my $tri_compo_kmer;
	if ($id_rates =~ m/(\w\w)\[(\w)\>\w\](\w\w)/xms){$tri_compo_kmer=join "", $1,$2,$3};
	my $value_of_compo = $composition_tri_penta{$tri_compo_kmer};
	
	#my $mut_tx = $count_value / ( $value_of_compo * $full_time_values);
	my $mut_tx;
	my $diviseur = ($value_of_compo / $number_tips) * $full_time_values;
	
	if ($count_value > 0 && $diviseur > 0 ){
		$mut_tx = $count_value / $diviseur;
	}else{
		$mut_tx = 0;
	};

	say {$out_compo_penta_tx} join  "\t", $id_rates , $mut_tx;
	$TxS_penta{$id_rates}=$mut_tx;
	
	$total_TXS_penta += $count_value;
};


close $out_compo_penta_tx;

	
	
my $out_matrix_compo_penta_txR = join "" , $outname,'_matrix_txR_penta.txt';
open my $out_compo_penta_txR, '>',  $out_matrix_compo_penta_txR;
say {$out_compo_penta_txR} join  "\t", "Mutation Types",$outname;


my $out_matrix_compo_penta_txE = join "" , $outname,'_matrix_txE_penta.txt';
open my $out_compo_penta_txE, '>',  $out_matrix_compo_penta_txE;
say {$out_compo_penta_txE} join  "\t", "Mutation Types",$outname;

## %TxS_penta

my %TxR_penta;
tie %TxR_penta, 'Tie::IxHash';
my %TxE_penta;
tie %TxE_penta, 'Tie::IxHash';

my $tot_TxR_penta_total=0;
my $tot_TxE_penta_total=0;

for my $id_TxS (keys %mix_matrix_penta){
	
	my $tri_compo_kmer;
	if ($id_TxS =~ m/(\w\w)\[(\w)\>\w\](\w\w)/xms){$tri_compo_kmer=join "", $1,$2,$3};
	my $value_of_compo = $composition_tri_penta{$tri_compo_kmer};
	
	
	
	my $val_TxS = $mix_matrix_penta{$id_TxS};
	my $id_TxR;
	if ($id_TxS =~ m/(\w+)\[(\w)\>(\w)\](\w+)/xms){$id_TxR=join "", $1,"[",$3,">",$2,"]",$4};
	my $val_TxR = $mix_matrix_penta{$id_TxR};
	
	my $true_revert;
	my $Rate_true_revert=0;
	if ($val_TxS > $val_TxR){$true_revert=$val_TxR}else{$true_revert=$val_TxS};
	
	my $diviseur = ($value_of_compo / $number_tips)* $full_time_values;
	
	if ($val_TxR > 0 && $diviseur > 0 ){
		$Rate_true_revert = $true_revert / $diviseur;
	};

	$TxR_penta{$id_TxS}=$Rate_true_revert;
	
	
	
	
	my $val_TxE=$val_TxS-$true_revert;	
	my $Rate_TxE=0;
	
	if ($val_TxE > 0 && $diviseur > 0 ){
		$Rate_TxE = $val_TxE / $diviseur;
	};

	$TxE_penta{$id_TxS}=$Rate_TxE;

	say {$out_compo_penta_txR} join  "\t", $id_TxS , $Rate_true_revert;
	say {$out_compo_penta_txE} join  "\t", $id_TxS , $Rate_TxE;
	
	$tot_TxR_penta_total+= $true_revert;
	$tot_TxE_penta_total+= $val_TxE;


}

close $out_compo_penta_txE;
close $out_compo_penta_txR;

my $TxS_penta_total = $total_TXS_penta/ ($full_time_values * $the_seq_length);
my $TxR_penta_total = $tot_TxR_penta_total/ ($full_time_values * $the_seq_length);
my $TxE_penta_total = $tot_TxE_penta_total/ ($full_time_values * $the_seq_length);



 
################################## Rate classique   ###################################################


my $TxS_classic = $N_total / ($full_time_values * $the_seq_length);
my $PING_TxS=$Ping_sum/($full_time_values * $the_seq_length);
my $Class_PING_TxE=$TxS_classic-$PING_TxS;

my $ratio_corrected = (sprintf  "%.3f", ($Ping_sum_select /$selected_count))*100;

my $out_report_Rates = join "" , $outname,'_Comparaison_sub_rate.txt';
open my $out_Ratesss, '>',  $out_report_Rates;
say {$out_Ratesss} join  "\t", "#ID","TxS_classic","TxS_mono","TxS_tri","TxS_penta","TxR_PING","TxR_mono","TxR_tri",
"TxR_penta","TxE_Class_PING","TxE_mono","TxE_tri","TxE_penta","Group","Genus","Family","Species","Prop_PING","Seq_length","Full_time",
"Time_tree","Ping_sum","Ping_sum_Norevertion","total_mut","Measuring_timescale","Selected_mut","Selected_Ping","Selected_ratio";

say {$out_Ratesss} join  "\t", 

	$outname,
	$TxS_classic,
	$TxS_mono,
	$TxS_total,
	$TxS_penta_total,
	$PING_TxS,
	$TxR_mono,
	$TxR_total,
	$TxR_penta_total,
	$Class_PING_TxE,
	$TxE_mono,
	$TxE_total,
	$TxE_penta_total,
	$group_rate,
	$genus_rate,
	$family_rate,
	$species_rate,
	$Ratio_PING,
	$the_seq_length,
	$full_time_values,
	$time_values,
	$Ping_sum,
	$Ping_sum_no_rest,
	$N_total,
	$Measuring_timescale,
	$selected_count,
	$Ping_sum_select,
	$ratio_corrected;
	
	

close $out_Ratesss;




######################################## Reformatage de matrix substitutions ########################

my $out_matrix_compo_txS_corrected = join "" , $outname,'_Corrected_matrix_txS_proportion.txt';
open my $out_compo_txS_corrected, '>',  $out_matrix_compo_txS_corrected;
say {$out_compo_txS_corrected} join  "\t", "Mutation Types",$outname;

my %TxS_corrected;
tie %TxS_corrected, 'Tie::IxHash';

## TxS_total:$TxS_total

for my $id_TxS (keys %TxS){
		
	my $val_TxS = $TxS{$id_TxS};
	my $val_TxS_CORRECTED =0;
	if ($TxS_total > 0 && $val_TxS > 0 ){
		$val_TxS_CORRECTED = ($val_TxS / $TxS_total ) * $N_total;
	}; 
	
	say {$out_compo_txS_corrected} join  "\t", $id_TxS , $val_TxS_CORRECTED;
	
}
close $out_compo_txS_corrected;





my $out_matrix_compo_txR_corrected = join "" , $outname,'_Corrected_matrix_txR_proportion.txt';
open my $out_compo_txR_corrected, '>',  $out_matrix_compo_txR_corrected;
say {$out_compo_txR_corrected} join  "\t", "Mutation Types",$outname;

my %TxR_corrected;
tie %TxR_corrected, 'Tie::IxHash';

for my $id_TxR (keys %TxR){
		
	my $val_TxR = $TxR{$id_TxR};
	my $val_TxR_CORRECTED =0;
	if ($TxR_total > 0 && $val_TxR > 0 ){
		$val_TxR_CORRECTED = ($val_TxR / $TxR_total ) * $N_total;
	}; 
	

	
	
	say {$out_compo_txR_corrected} join  "\t", $id_TxR , $val_TxR_CORRECTED;
	
}
close $out_compo_txR_corrected;




my $out_matrix_compo_txE_corrected = join "" , $outname,'_Corrected_matrix_txE_proportion.txt';
open my $out_compo_txE_corrected, '>',  $out_matrix_compo_txE_corrected;
say {$out_compo_txE_corrected} join  "\t", "Mutation Types",$outname;

my %TxE_corrected;
tie %TxE_corrected, 'Tie::IxHash';

for my $id_TxE (keys %TxE){
		
	my $val_TxE = $TxE{$id_TxE};
	
	my $val_TxE_CORRECTED =0;
	
	if ($TxE_total > 0 && $val_TxE > 0 ){
		$val_TxE_CORRECTED = ($val_TxE / $TxE_total ) * $N_total;
	}; 
		
	say {$out_compo_txE_corrected} join  "\t", $id_TxE , $val_TxE_CORRECTED;
	
}

close $out_compo_txE_corrected;






my $out_matrix_compo_penta_txS_corrected = join "" , $outname,'_Corrected_matrix_txS_penta_proportion.txt';
open my $out_compo_penta_txS_corrected, '>',  $out_matrix_compo_penta_txS_corrected;
say {$out_compo_penta_txS_corrected} join  "\t", "Mutation Types",$outname;

my %TxS_penta_corrected;
tie %TxS_penta_corrected, 'Tie::IxHash';

for my $id_TxS_penta (keys %TxS_penta){
		
	my $val_TxS_penta = $TxS_penta{$id_TxS_penta};
	my $val_TxS_penta_CORRECTED = 0; 
	
	if ($TxS_penta_total > 0 && $val_TxS_penta > 0 ){
		$val_TxS_penta_CORRECTED = ($val_TxS_penta / $TxS_penta_total ) * $N_total;
	}; 
	
	say {$out_compo_penta_txS_corrected} join  "\t", $id_TxS_penta , $val_TxS_penta_CORRECTED;
	
}

close $out_compo_penta_txS_corrected;



my $out_matrix_compo_penta_txR_corrected = join "" , $outname,'_Corrected_matrix_txR_penta_proportion.txt';
open my $out_compo_penta_txR_corrected, '>',  $out_matrix_compo_penta_txR_corrected;
say {$out_compo_penta_txR_corrected} join  "\t", "Mutation Types",$outname;

my %TxR_penta_corrected;
tie %TxR_penta_corrected, 'Tie::IxHash';

for my $id_TxR_penta (keys %TxR_penta){
		
	my $val_TxR_penta = $TxR_penta{$id_TxR_penta};
	my $val_TxR_penta_CORRECTED = 0; 
	
	if ($TxR_penta_total > 0 && $val_TxR_penta > 0 ){
		$val_TxR_penta_CORRECTED = ($val_TxR_penta / $TxR_penta_total ) * $N_total;
	}; 
	
	say {$out_compo_penta_txR_corrected} join  "\t", $id_TxR_penta , $val_TxR_penta_CORRECTED;
	
}

close $out_compo_penta_txR_corrected;



my $out_matrix_compo_penta_txE_corrected = join "" , $outname,'_Corrected_matrix_txE_penta_proportion.txt';
open my $out_compo_penta_txE_corrected, '>',  $out_matrix_compo_penta_txE_corrected;
say {$out_compo_penta_txE_corrected} join  "\t", "Mutation Types",$outname;

my %TxE_penta_corrected;
tie %TxE_penta_corrected, 'Tie::IxHash';

for my $id_TxE_penta (keys %TxE_penta){
		
	my $val_TxE_penta = $TxE_penta{$id_TxE_penta};
	my $val_TxE_penta_CORRECTED = 0; 
	
	if ($TxE_penta_total > 0 && $val_TxE_penta > 0 ){
		$val_TxE_penta_CORRECTED = ($val_TxE_penta / $TxE_penta_total ) * $N_total;
	}; 
	
	say {$out_compo_penta_txE_corrected} join  "\t", $id_TxE_penta , $val_TxE_penta_CORRECTED;
	
}


close $out_compo_penta_txE_corrected;


##########################################   Matrix normalization #####################################


my $outname_final_normalized = join "" , $outname,'_matrix_normalized_line.txt';
open my $out_fin_normalized, '>',  $outname_final_normalized;
say {$out_fin_normalized} join  "\t", "Mutation Types",$outname;

my $outname_final_normalized_positiv = join "" , $outname,'_matrix_normalized_positiv_line.txt';
open my $out_fin_normalized_positiv, '>',  $outname_final_normalized_positiv;
say {$out_fin_normalized_positiv} join  "\t", "Mutation Types",$outname;

my $outname_final_normalized_negativ = join "" , $outname,'_matrix_normalized_negativ_line.txt';
open my $out_fin_normalized_negativ, '>',  $outname_final_normalized_negativ;
say {$out_fin_normalized_negativ} join  "\t", "Mutation Types",$outname;

my $test_tot =0;
my $test_tot_norm =0;
my $full_prop =0;
my $full_prop2 =0;
my $iii=0;



for my $id_rates (keys %mix_matrix){
	my $tri_compo_kmer;
	if ($id_rates =~ m/(\w)\[(\w)\>\w\](\w)/xms){$tri_compo_kmer=join "", $1,$2,$3};
	my $value_of_compo = mean(@{$composition_tri{$tri_compo_kmer}});
	
	
	my $prop = ($value_of_compo/ $full_tri_compo); 
	
	$full_prop +=  $prop ;

}





for my $id_rates (keys %mix_matrix){
	my $count_value = $mix_matrix{$id_rates};
	$test_tot  += $count_value;
	my $tri_compo_kmer;
	if ($id_rates =~ m/(\w)\[(\w)\>\w\](\w)/xms){$tri_compo_kmer=join "", $1,$2,$3};
	my $value_of_compo = mean(@{$composition_tri{$tri_compo_kmer}});
	
	## $tri_compo_kmer
	## $full_tri_compo
	## $value_of_compo
	## $N_total
	
	my $normalizer;
	
	my $prop = ($value_of_compo/ $full_tri_compo)/$full_prop; 
	
	## $prop
	$iii++;
	
	$full_prop2 +=  $prop ;
	if ($value_of_compo > 0 && $full_tri_compo > 0 ){
		$normalizer  = (($value_of_compo/ $full_tri_compo)/$full_prop )* $N_total;
	}else{
		$normalizer = 0;
	};
	
	my $count_value_norm = $count_value - $normalizer;
	
	
	$test_tot_norm  += $normalizer;
	##  $count_value
	##  $normalizer
	##  $count_value_norm
	
	$count_value_norm = sprintf  "%.0f", $count_value_norm;
	
	say {$out_fin_normalized} join  "\t", $id_rates , $count_value_norm;
	
	if ($count_value_norm > 0 ){say {$out_fin_normalized_positiv} join  "\t", $id_rates , $count_value_norm;}
	else {say {$out_fin_normalized_positiv} join  "\t", $id_rates , 0;}
	
	if ($count_value_norm < 0 ){say {$out_fin_normalized_negativ} join  "\t", $id_rates , $count_value_norm;}
	else {say {$out_fin_normalized_negativ} join  "\t", $id_rates , 0;}
	
	
}

## $iii
## $full_prop
## $full_prop2

## $test_tot_norm 
## $N_total
## $test_tot


close $out_fin_normalized;		
close $out_fin_normalized_positiv;	
close $out_fin_normalized_negativ;	


my $outname_final_pinta_normalized = join "" , $outname,'_matrix_pinta_normalized_3086_line.txt';
open my $out_fin_pinta_normalized, '>',  $outname_final_pinta_normalized;
say {$out_fin_pinta_normalized} join  "\t", "Mutation Types",$outname;

my $outname_final_pinta_normalized_positiv = join "" , $outname,'_matrix_pinta_normalized_3086_positiv_line.txt';
open my $out_fin_pinta_normalized_positiv, '>',  $outname_final_pinta_normalized_positiv;
say {$out_fin_pinta_normalized_positiv} join  "\t", "Mutation Types",$outname;

my $outname_final_pinta_normalized_negativ = join "" , $outname,'_matrix_pinta_normalized_3086_negativ_line.txt';
open my $out_fin_pinta_normalized_negativ, '>',  $outname_final_pinta_normalized_negativ;
say {$out_fin_pinta_normalized_negativ} join  "\t", "Mutation Types",$outname;


for my $id_rates (keys %matrix_penta){
	my $count_value = $matrix_penta{$id_rates};
	
	my $tri_compo_kmer;
	if ($id_rates =~ m/(\w\w)\[(\w)\>\w\](\w\w)/xms){$tri_compo_kmer=join "", $1,$2,$3};
	my $value_of_compo = $composition_tri_penta{$tri_compo_kmer};
	
	my $normalizer;
	
	if ($value_of_compo > 0 && $full_tri_compo > 0 ){
		$normalizer  = ($value_of_compo / $full_penta_compo ) * $N_total;
	}else{
		$normalizer = 0;
	};
	
	my $count_value_norm = $count_value - $normalizer;
	
	say {$out_fin_pinta_normalized} join  "\t", $id_rates , $count_value_norm;
	
	if ($count_value_norm > 0 ){say {$out_fin_pinta_normalized_positiv} join  "\t", $id_rates , $count_value_norm;}
	else {say {$out_fin_pinta_normalized_positiv} join  "\t", $id_rates , 0;}
	
	if ($count_value_norm < 0 ){say {$out_fin_pinta_normalized_negativ} join  "\t", $id_rates , $count_value_norm;}
	else {say {$out_fin_pinta_normalized_negativ} join  "\t", $id_rates , 0;}
	
	
}

close $out_fin_pinta_normalized;	
close $out_fin_pinta_normalized_positiv;	
close $out_fin_pinta_normalized_negativ;		
	



# ########################################    Functions    ############################################


sub mean {
  return sum(@_)/@_;
}

sub SD{
	my $set = shift ;
	my @set= @$set;

	my $m = mean(@set);
	my $n = scalar @set;

	my $sum; 

	foreach my $i (@set){
	
		my $partial = $i - $m;
		$partial *= $partial;
		$sum += $partial;

	}	
	
	my $sd = "Na";
	
	if ($n >1 ){
	
		my $int_sd = $sum / ($n-1);
		$sd = sqrt($int_sd);
	
	};
	
	return $sd;

};

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


sub read_ID {

	my $infile = shift ;	

	## Reading input file: $infile
	## Elapsed time |===[%]

	open my $in, '<', $infile;

	my %ID_refs;
	my %ID_ref_rev;

	LINE:
		
	while (my $line = <$in> ){
		chomp $line;
		
		if ($line =~ m/(\d+)\ ([^\|]+\|[^\,]+)/xms) {
			$ID_refs{$1}=$2;
			$ID_ref_rev{$2}=$1;
		}
		
		next LINE;
		
	};
	
	return (\%ID_refs,\%ID_ref_rev);
};




sub read_length {

	my $infile = shift ;	

	## Reading input file: $infile
	## Elapsed time |===[%]

	open my $in, '<', $infile;

	my %ID_len;

	LINE:
		
	while (my $line = <$in> ){
		chomp $line;
		
		if ($line =~ m/([^\t]+)\t(\d+)/xms) {
			$ID_len{$1}=$2;
		}
		
		next LINE;
		
	};
	
	return %ID_len;
};





sub read_tree{

	## Generate tree hash

	my $tree_ref_file = shift;
		
	open my $in, '<', $tree_file;		  
	my $tree_ref = <$in>;
	chomp $tree_ref;
	$tree_ref =~ s/.*\[&R\]//;
	my $treeT  = $tree_ref;
	$treeT =~ s/[^\(\)\,]+//g;
	$treeT =~ s/\),\(//g;
	$treeT =~ s/,//g;
	$treeT =~ s/\)//g;
	my $Ngen = (length $treeT) + 2;

	my @Life_rate;
	my %IDfeuilles;
	my %tree_hash;
	my %result_hash;
	my $e0 = 0;
	my $id01 = "N0A0";


	my $TOTAL_A_count = 0;
	my $TOTAL_T_count = 0;
	my $TOTAL_C_count = 0;
	my $TOTAL_G_count = 0;
	my $pop =0;



	## $tree_ref

	if ($tree_ref =~ m/.+\[\&.+\=\"(.+)\"\];/xms){


							
		my $date=0;
		my $veriRate =0;
		my $veri_seq = $1;
			
		$tree_ref =~ s/(.+)\[\&.+\=\"(.+)\"\];/$1/g;

		##  $veri_seq 

		$tree_hash{$e0}{$id01}= {
			"tree" => $tree_ref,
			"Value"=> "0",
			"seq"=> $veri_seq,
			"Nature" => "Last_Ancestor",
			"Date"=> $date,
			"Rate"=> $veriRate
		
		};
							
		$result_hash{$id01}= {
			"tree" => $tree_ref,
			"Value"=> "0",
			"seq"=> $veri_seq,
			#"anc_seq" => $tree_hash{"1"}{""}{"seq"},
			"Nature" => "Last_Ancestor",
			"Date"=> $date,
			"Rate"=> $veriRate
		
		};

							
	};
	
	
	## $tree_ref	
	### $Ngen	

	SPLIT:
	for (my $y = 0 ; $y < $Ngen; $y +=1 ){
		## $y
		IDS:
		for my $id (keys %{$tree_hash{$y}}){
			## $id
			my $tree = $tree_hash{$y}{$id}{"tree"};

			my $tree1= $tree;
			my $PATTERN = $RE{balanced}{-parens=>'()'};
			my @matches = ();
			
			## $tree1		
			INPUT:
			while (length $tree) {
				$tree =~ s{\A [^\w\(]+ }{}x;  ## skip leading junk
				if ($tree =~ s{ ($PATTERN\[[^\]]+\]\:\[[^\]]+\][^\,^\)]+) }{}x ) {
					push @matches, $1;
						
					next INPUT;
				}
					
				last INPUT;
			
			};

			my @values;
				
				if ( $tree =~ m/\,/xms ) {@values = split (/\,/, $tree);
			
					}
				else { push @values, $tree};
			
			## @matches
			## @values

			my $i =0;
			my $e = $y +1;

				
			while ( @matches ){
				my $ext = shift @matches;
				my ($bout_name) = $id =~ m/N\d+(A\d+)/xms;
				my $new_id = join "","N$e",$bout_name,$i; 

				if($ext =~ m/\((.+)\)\[&[^\=]+\=\"(\w+)\"\]\:\[.rate\=([^\]]+)\]([^\,^\)]+)/xms){
					
					my $date=lc($4);
						
					my $veri_seq = $2;
					my $veriTree =$1;
					my $veriRate =lc ($3);
						
					if ( $date =~ m/(\d+\.?\d+?)e([\-\+])*(\d*)/xms ){

					my $numprem = $1;
						my $puisss = $3;
						my $tutu = $2;
						
						if (defined $tutu){
							if ($2 =~ m/\-/xms ) {
								$veriRate= $date / (10 ** $puisss);
							}else {
								$veriRate= $date * (10 ** $puisss);
							};
						};
					};
						
											
					if ( $veriRate =~ m/(\d+\.?\d+?)e([\-\+])*(\d*)/xms ){
						
						
						my $numprem = $1;
						my $puisss = $3;
						my $tutu = $2;
						
						if (defined $tutu){
							if ($2 =~ m/\-/xms ) {
								$veriRate= $numprem / (10 ** $puisss);
							}else {
								$veriRate= $numprem * (10 ** $puisss);
							};
						};
					};
						
						my $In_mother_seq = $tree_hash{$y}{$id}{"seq"};
						
						
						my $ixt= $ext;
															
						$ixt =~ s/\(.+\)(\[&[^\=]+\=\"\w+\"\]\:\[.rate\=[^\]]+\][^\,^\)]+)/$1/g;
						
						my $axt= $ext;
					 
						$axt =~ s/\(.+\)(\[&[^\=]+\=\")\w+(\"\]\:\[.rate\=[^\]]+\][^\,^\)]+)/$1$new_id$2/g;
						
						if ( $id eq $id01 ){ $In_mother_seq =$veri_seq };
						
						my $mask= $veri_seq ^ $In_mother_seq;
						
						my @n_mut = ($mask =~ /[^\0]/g);
						
						my $n_mut = @n_mut;
						


						$tree_hash{$e}{$new_id}= {
							"tree" => $veriTree,
							"seq"=> $veri_seq,
							"Nature"=>"Ancestor",
							"Date"=> $date,
							"Rate"=> $veriRate,
							"Nmute"=>$n_mut
						};
						
						$result_hash{$new_id}= {
							"anc_seq" => $tree_hash{$y}{$id}{"seq"},
							"seq"=> $veri_seq,
							"Nature"=>"Ancestor",
							"Date"=> $date,
							"Rate"=> $veriRate,
							"Nmute"=>$n_mut
						};

						my $A_count = $veri_seq =~ tr/A//;
						my $T_count = $veri_seq =~ tr/T//;
						my $C_count = $veri_seq =~ tr/C//;
						my $G_count = $veri_seq =~ tr/G//;
						
						$TOTAL_A_count = $TOTAL_A_count + $A_count;
						$TOTAL_T_count = $TOTAL_T_count + $T_count ;
						$TOTAL_C_count = $TOTAL_C_count + $C_count;
						$TOTAL_G_count = $TOTAL_G_count + $G_count;
							$i++;
						};

				};
			
				
			
			my @values2=@values;
			while ( @values ){
				
				## 2
				## $id
				my $ext = shift @values;
				my ($bout_name) = $id =~ m/N\d+(A\d+)/xms;
				my $new_id = join "","N$e",$bout_name,$i; 
				
				## $ext
	#			if ($ext =~ m/\)\[\&.+\=\"(.+)\"\];/xms){
					
	#				my $date=0;
	#								my $veriRate =0;
					
	#				$tree_hash{$e0}{$new_id}= {
	#					"seq"=> $1,
	#					"Nature" => "Last_Ancestor",
	#					"Date"=> $date,
	#					"Rate"=> $veriRate
	#					};
					
									
	#				$result_hash{$new_id}= {
	#					"seq"=> $1,
	#					"anc_seq" => $tree_hash{$y}{$id}{"seq"},
	#					"Nature" => "Last_Ancestor",
	#					"Date"=> $date,
	#					"Rate"=> $veriRate
	#				};
						
	#				$i++;
						
			#	}

				
				if($ext =~ m/(\d+)(\[&[^\=]+\=\"(\w+)\"\]\:\[.rate\=([^\]]+)\])([^\,^\)]+)/xms){
				
					my $date=lc($5);
					my $identifier = $2;
					my $veri_seq = $3;
					my $numerous = $1;
					my $veriRate =lc ($4);
					
					if ( $date =~ m/(\d+\.?\d+?)e([\-\+])*(\d*)/xms ){
						my $numprem = $1;
						my $puisss = $3;
						my $tutu = $2;
						
						if (defined $tutu){
							if ($2 =~ m/\-/xms ) {
								$date= $date / (10 ** $puisss);
							}else {
								$date= $date * (10 ** $puisss);
							};
						};
					};
					
										
					if ( $veriRate =~ m/(\d+\.?\d+?)e([\-\+])*(\d*)/xms ){
						## $date
						## lala
						my $numprem = $1;
						my $puisss = $3;
						my $tutu = $2;
						
						if (defined $tutu){
							if ($2 =~ m/\-/xms ) {
								$veriRate= $numprem / (10 ** $puisss);
							}else {
								$veriRate= $numprem * (10 ** $puisss);
							};
						};
					};
						##3
								## $y
					## $id
					## $new_id

					
					$tree_hash{$e}{$new_id}= {
						"tree" => "",
						"seq"=> $veri_seq,
						"Nature" => $ID_ref{$numerous},
						"Date"=> $date,
						"Rate"=> $veriRate
					};
					
					
					$result_hash{$new_id}= {
						"anc_seq" => $tree_hash{$y}{$id}{"seq"},
						"seq"=> $veri_seq,
						"Nature" => $ID_ref{$numerous},
						"Date"=> $date,
						"Rate"=> $veriRate,
						"IDENT"=> $identifier
					};
					
					my $A_count = $veri_seq =~ tr/A//;
					my $T_count = $veri_seq =~ tr/T//;
					my $C_count = $veri_seq =~ tr/C//;
					my $G_count = $veri_seq =~ tr/G//;
					
					$TOTAL_A_count = $TOTAL_A_count + $A_count;
					$TOTAL_T_count = $TOTAL_T_count + $T_count ;
					$TOTAL_C_count = $TOTAL_C_count + $C_count;
					$TOTAL_G_count = $TOTAL_G_count + $G_count;
					
					
					push @Life_rate, $veriRate;

					$IDfeuilles{$new_id}=$ID_ref{$numerous};
						
					$i++;
						
				};


			};
			
					
		};
		
	};


	return ( \%tree_hash, \%result_hash, \%IDfeuilles);

};




sub read_info {

	my %Read_info;
	my $infile = shift ;	

	open my $in, '<', $infile;

	LINE:
		
	while (my $line = <$in> ){
		
		chomp $line;
		my @infos=split (/\_/, $line);
		## @infos
		my $Species = lc (shift @infos);
		## $Species
		## @infos
		
		$Species=~s/human\ //g;
		$Species =~ s/\ /\_/g;
		$Species =~ s/\-//g;
		## $Species
		
		my $Genus;
		my $family;
		my $strain;
		my $group ;
		
		
		## y:length(@infos)
		## x:scalar(@infos)
		if (scalar(@infos) == 4 ){
		$Genus = shift  @infos;
		$family = shift @infos;
		$strain = shift @infos;
		$group = shift @infos;
	}else{
		$Genus = shift  @infos;
		$family = shift @infos;
		$strain = "Na";
		$group = shift @infos;
		
		};
		$Read_info{$Species}={
		
			"Group" => $group,
			"Genus" => $Genus,
			"Family" => $family,
			"Strain" => $strain,
		
		};
	
	
	}		
	
	return %Read_info ;

};

