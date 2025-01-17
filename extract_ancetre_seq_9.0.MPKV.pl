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

my $life_seq = shift;
my $outname_short = shift;


##########
my @exclusion;
my $exclusion_file = shift;
if( defined $exclusion_file ){
open my $in_exclusion, '<', $exclusion_file;



LINE:
		
while (my $line_exclusion = <$in_exclusion> ){

	chomp $line_exclusion;

	push @exclusion, $line_exclusion;

};
}

#########





my ($outname) = $tree_file =~ m/(.*)\.\w+/g;


my ($filespecies) = $tree_file =~ m/(.*)\.fa.*\.tree.*/g; 


my %Ancestor_hash = read_Anc( $ID_file );



read_Anc_life( $life_seq );


## x:keys(%Ancestor_hash)

## %ID_ref
my %Info_hash = read_info( $INFO_table);

my @keys= keys %Ancestor_hash;
### @keys

my ($tree_hash, $result_hash, $ID_feuilles) = read_tree( $tree_file);
my %tree_hash = %$tree_hash;
my %result_hash = %$result_hash;
my %ID_feuilles = %$ID_feuilles; 


my %Info_length = read_length( $INFO_size);
### mpx9

### $filespecies
### $outname
### $outname_short

my $the_seq_length;


my @tmp_length = keys %Info_length;

if (not defined $the_seq_length){

	my $tarttuf = $outname_short;
	$tarttuf =~ s/\_iter\-N\d+.+//g;
	my @tmp_grepage = grep (/${tarttuf}/, @tmp_length);
	my $tmp_grepage = $tmp_grepage[0];

	if (defined $tmp_grepage){
	
		$the_seq_length = $Info_length{$tmp_grepage};
	
	}else{
		
		$tarttuf =~ s/\w\w\-iT\_//g;
		my @speciespecies = split "_", $tarttuf;
		
		my $tarttuf2 = "TATA";
		
		LENELENELEN:
		while (not defined $the_seq_length){
			
			$tarttuf =~ s/\_[^\_]*?$//;
			
			
			my @tmp_grepage_recherch = grep (/${tarttuf}/, @tmp_length);
			
			my $tmp_grepage_recherch = $tmp_grepage_recherch[0];

			if (defined $tmp_grepage_recherch){
	
				$the_seq_length = $Info_length{$tmp_grepage_recherch};
	
			};
			 
			if ($tarttuf eq $tarttuf2){ last LENELENELEN };
			$tarttuf2 = $tarttuf;
			
			next LENELENELEN;
		
		};
		
		
		if (not defined $the_seq_length){$the_seq_length = length ($result_hash{"N1A00"}{"seq"})};
		
	};

	
};




my %TOtal_mutation;
my $N_total = 0;
my $full_time_values=1;
my $kj=0;		
my $group = "dsDNA";
my $genus= "Poxviridae";
my $family= "Orthopoxvirus";
my $species = "MPXV";					
					
my $group_rate =  $group ;
my $genus_rate = $genus;
my $family_rate = $family;
my $species_rate = $species;	



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
	my $seq = $result_hash{$id}{"seq"};
	my $nat = $result_hash{$id}{"Nature"};
	if ( grep /$nat/ , @exclusion){next TAST};
	### $id
	
	if( $id eq "N0A0" ){next TAST};
	my %res;
	my $time_data = 0;
	
	## $id

	
	### $nat

	
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
my $time_data = 0;

## $seq_rf
## @DNA_ref

my $total_mut=0;



TURN:
while(@DNA_ref){
	my $nec = shift @DNA_ref;
	## $nec
	## $i
	
	if ($nec =~ m/\w/xms){ 
	for my $id (keys %{$res{$i}}){
				##$id
		
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
				my $down1;
				
				
				if ( $u1 < length $seq_rf ){
					
					$down1 = substr ($seq_rf, $d1, 1);	
					DOWN1:
					while ($down1 =~ m/\-/xms ){
						$d1++;
						$down1 = substr ($seq_rf, $d1, 1);
						
					};
				
				}else{
				
					$down1 = "A";
				
				};
				
				my $d2 = $d1+1;
				my $down2;
				
				
				if ( $d2 < length $seq_rf ){
					
					$down2 = substr ($seq_rf, $d2, 1);
				
					DOWN2:
					while ($down2 =~ m/\-/xms ){
						$d2++;
						$down2 = substr ($seq_rf, $d2, 1);
					};
				
				}else{ 
					$down2 = "A";
					
				};
						 				

				
				## $group
				
				
				$N_total++;
				## $id
				
				say {$out} join "\t",  $id,$nec, $i+2, $res{$i}{$id}{"Nucl"}, $res{$i}{$id}{"Codon"}, (join "",$up2,$up1,"[",$nec,">",$res{$i}{$id}{"Nucl"},"]",$down1,$down2),
				(join "",$nec,">",$res{$i}{$id}{"Nucl"}),$statut,$time_data,$full_time_values,0,$realseqlength, $group,$genus,$family,$species,$filespecies;
				
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
		say {$out_fin_single_compil} join "\t",$outname,$tres,$matrix_of_singl{$tres},$group,$genus,$family,$species,$the_seq_length,$full_time_values,0;
		}
	else {say {$out_fin_single_compil} join "\t", $outname,$tres, 0,$group,$genus,$family,$species,$the_seq_length,$full_time_values,0;
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












# ################## Analyse Sub Rev Evo rates mono context ###############################


# ################ matrix compo 

### test 0


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



### test 1

	
my $k = 'A,T,G,C';
my $kmero_compo_mono;

my %composition_tri_mono;
tie %composition_tri_mono, 'Tie::IxHash';



for  my $kmer_compo_mono (glob"{$k}"){
	$composition_tri_mono{$kmer_compo_mono}=0;
}

my $full_mono_compo  = 0;

### test 2


for my $id_compo_mono (keys %seq_for_compo_mono){

	$seq_compo_mono = $seq_for_compo_mono{$id_compo_mono};
	
	for my $kmer_compo_mono (keys %composition_tri_mono){
		
		my $countage_compo_mono = () = $seq_compo_mono =~ /$kmer_compo_mono/g;

		$composition_tri_mono{$kmer_compo_mono}+= $countage_compo_mono;
		$full_mono_compo +=  $countage_compo_mono;
		
	}
}

### test 3

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
	
	say {$out_fin_single_compil_txs} join "\t",$outname,$id_TxS_mono,$val_TxS,$group,$genus,$family,$species,$the_seq_length,$full_time_values,0;

};
	


close $out_fin_single_compil_txs;












my $out_matrix_compo_txR_mono = join "" , $outname,'_matrix_txR_mono.txt';
open my $out_compo_txR_mono, '>',  $out_matrix_compo_txR_mono;
say {$out_compo_txR_mono} join  "\t", "Mutation Types",$outname;


my $out_matrix_compo_txE_mono = join "" , $outname,'_matrix_txE_mono.txt';
open my $out_compo_txE_mono, '>',  $out_matrix_compo_txE_mono;
say {$out_compo_txE_mono} join  "\t", "Mutation Types",$outname;



					
					


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
	0,
	$PING_TxS,
	$TxR_mono,
	$TxR_total,
	0,
	$Class_PING_TxE,
	$TxE_mono,
	$TxE_total,
	0,
	$group_rate,
	$genus_rate,
	$family_rate,
	$species_rate,
	$Ratio_PING,
	$the_seq_length,
	$full_time_values,
	0,
	$Ping_sum,
	$Ping_sum_no_rest,
	$N_total,
	1,
	$selected_count,
	$Ping_sum_select,
	$ratio_corrected;
	
	

close $out_Ratesss;




	



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
		
	open my $in, '<', $tree_ref_file;	
	my $tree_ref;
	
	LINE:
	while (my $tree_line = <$in> ){
		
		chomp $tree_line;
		if ($tree_line =~ m/tree\ tree\_1\ \=\ \[\&R\]\ +(.+)/xms) {
			$tree_ref=$1;
		}
		
		next LINE;
	
	}
	
	
	$tree_ref =~ s/(Node\d+)\/(\d+)/$1\@$2/g;
	$tree_ref =~ s/\///g;
	$tree_ref =~ s/@/\//g;
	$tree_ref =~ s/\-//g;
	$tree_ref =~ s/\'//g;
	
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


	
	### $tree_ref


	
	if($tree_ref =~ m/\((.+\)([^\)]*)\:[^\,]+?),[^\,]+\;/xms){
		
		### juju
		
		my $tree_ref = $1;
		my $ident=$2;
		my $veri_seq="NA";
		
		### $ident
		
			my $five_first 	=$veri_seq;
			if ($five_first ){$five_first =~ s/(^.{5}).*$/$1/;};
				
				### $five_first
		
		if ($ident){
			
			if ( $ident =~ m/(.+)\/.+/xms){
						
				$ident=$1;
		
			};
		
			$veri_seq  = $Ancestor_hash{$ident};

		};
		
		$tree_hash{$e0}{$id01}= {
			"tree" => $tree_ref,
			"seq"=> $veri_seq,
			"Nature" => "Last_Ancestor"
		
		};
							
		$result_hash{$id01}= {
			"tree" => $tree_ref,
			"seq"=> $veri_seq,
			"Nature" => "Last_Ancestor"
		
		};



	};
		

	SPLIT:
	for (my $y = 0 ; $y < $Ngen; $y +=1 ){
		
		## $y
		
		IDS:
		for my $id (keys %{$tree_hash{$y}}){
			
			### $id

			my $tree = $tree_hash{$y}{$id}{"tree"};
			
		
			
			my $tree1= $tree;
			my $PATTERN = $RE{balanced}{-parens=>'()'};
			my @matches = ();
			
			INPUT:
			while (length $tree) {
				#$tree =~ s{\A [^\w\(]+ }{}x;  ## skip leading junk
				if ($tree =~ s{ ($PATTERN[^\:\)]*\:[^\,^\)]+) }{}x ) {
					push @matches, $1;
						
					next INPUT;
				}
					
				last INPUT;
			
			};

			my @values;
				
				if ( $tree =~ m/\,/xms ) {@values = split (/\,/, $tree);
			
					}
				else { push @values, $tree};
			
			my $i =0;
			my $e = $y +1;
			
			
	
			
				
			while ( @matches ){
				my $ext = shift @matches;
				## $ext
				
				my ($bout_name) = $id =~ m/N\d+(A\d+)/xms;
				my $new_id = join "","N$e",$bout_name,$i; 

				if($ext =~ m/\((.+)\)\[\&label\=\"(Node[\d\/]+)\"\]\:([\d\.E\-]+)/xms){
					
					my $veriTree =$1;
					my $ident=$2;
				
					if ( $ident =~ m/(.+)\/.+/xms){
						
						$ident=$1;
					};

					my $veri_seq  = $Ancestor_hash{$ident};
				
						## $ext 
							### $new_id
							## $veriTree
							### $ident
							## Acestor
							
							###X:length $veri_seq
							
					
				my $five_first = $veri_seq;
				if ($five_first ){$five_first =~ s/(^.{5}).*$/$1/;};
				
				### $five_first
						
					$tree_hash{$e}{$new_id}= {
						"tree" => $veriTree,
						"seq"=> $veri_seq,
						"Nature"=>"Ancestor",
					};
						
					$result_hash{$new_id}= {
						"anc_seq" => $tree_hash{$y}{$id}{"seq"},
						"seq"=> $veri_seq,
						"Nature"=>"Ancestor",

					};

					$i++;
				};

			};
			
			
			my @values2=@values;
			while ( @values ){
				
				
				
				my $ext = shift @values;
				
				## $ext
				
				my ($bout_name) = $id =~ m/N\d+(A\d+)/xms;
				my $new_id = join "","N$e",$bout_name,$i; 
				
				if($ext =~ m/([^\(]+)\:[\d+\.E\-]+/xms){
				
					my $identifier = $1;
					$identifier =~ s/\|/\_/g;
					$identifier =~ s/\//\_/g;
			
					my $veri_seq  = $Ancestor_hash{$identifier};
					
							## $ext 
							### $new_id
							### $identifier
							###X:length($veri_seq)
							
									my $five_first = $veri_seq;
					if ($five_first ){$five_first =~ s/(^.{5}).*$/$1/;};
				
				### $five_first

					
					$tree_hash{$e}{$new_id}= {
						"tree" => "",
						"seq"=> $veri_seq,
						"Nature" => $identifier
					};
					
					
					$result_hash{$new_id}= {
						"anc_seq" => $tree_hash{$y}{$id}{"seq"},
						"seq"=> $veri_seq,
						"Nature" => $identifier,
						"IDENT"=> $identifier
					};
					

					$IDfeuilles{$new_id}=$identifier;
						
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



sub read_Anc {

	my %seqset;
	my $infile = shift ;
	open my $in, '<', $infile;

	## START
	LONE:		
	while (my $line = <$in> ){
		chomp $line;
		 ## $line
		if ($line =~ m/Node\d+/xms){
			$line =~ s/\|/\_/g;
			
			$line =~ s/\//\_/g;
			my @infox = split ( "\t",$line);
			
			my $ID = $infox[0];
			my $Pose = $infox[1];
			my $seq =$infox[2];
			
			$seqset{$ID} .= $seq;

		};


	};

	return %seqset ;

};



sub read_Anc_life {

my $infile = shift ;	

	open my $in, '<', $infile;

	my $seq_id;
	my $seq;

	LINE:
		
		while (my $line = <$in> ){
		chomp $line;

		# at each '>' char...
		if ($line =~ m/>(.+)/xms) {
			
		 # add current seq to hash  (if any)
		 if($seq) {
			$Ancestor_hash{$seq_id} = $seq;
			$seq = q{};
			}
		
			# extract new seq_id
			$seq_id = $1;
					$seq_id =~ s/\///g;
		$seq_id =~ s/\-//g;
			$seq_id =~ s/\|/\_/g;
			$seq_id =~ s/\//\_/g;
			
			
			next LINE;
			
		}
		# elongate current seq (seq can be broken sevetal lines)
		$seq .= $line;
	}

	#add last seq to hash (if any)
	if($seq) {
		
		$Ancestor_hash{$seq_id} = $seq;
		$seq = q{};
	};

	close $in;

};
