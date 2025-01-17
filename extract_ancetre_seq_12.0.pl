#!/usr/bin/env perl
use warnings;
use strict;
use Modern::Perl '2011';
use autodie;
use Smart::Comments '###';
use List::AllUtils 'mesh';
use Algorithm::Numerical::Shuffle qw/shuffle/;
use List::Util qw(sum);
use List::Util qw(max);
use List::Util qw(min);
use Tie::IxHash;
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

use Bag::Similarity::Cosine;
#use Statistics::Basic qw(:all nofill);

my $tree_file = shift;
my $ID_file = shift;
my $INFO_table= shift;
my $INFO_size= shift;




my ($outname) = $tree_file =~ m/(.*)\.\w+/g;
my ($outname_short) = $tree_file =~ m/Uniq\_true\_tree\_(.+)\_sub_selected\_iteration\_\d+\.tree/g;
my $iteratio = 0;
if ($tree_file =~ m/iteration_(\d+)/g){$iteratio = $1};


my ($ID_ref,$Id_ref_rev) = read_ID( $ID_file );
my %Id_ref_rev = %$Id_ref_rev;
my %ID_ref = %$ID_ref;

## %ID_ref
my %Info_hash = read_info( $INFO_table);
my ($tree_hash, $result_hash, $ID_feuilles) = read_tree( $tree_file);
my %tree_hash = %$tree_hash;
my %result_hash;
tie %result_hash, 'Tie::IxHash';
%result_hash = %$result_hash;
my %ID_feuilles = %$ID_feuilles; 


my %Info_length = read_length( $INFO_size);


### $tree_file
### $outname
### $outname_short

my $the_seq_length = $Info_length{$outname_short} ;


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
my %ID_times;

#my $filespecies="NA";
$Info_hash{"NA"}{"Group"}="Ambigous";
$Info_hash{"NA"}{"Genus"}="Ambigous";
$Info_hash{"NA"}{"Family"}="Ambigous";
$Info_hash{"NA"}{"Strain"}="Ambigous";
	
	

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
	$full_time_values = $full_time_values + $result_hash{$id}{"Date"} ;#if $id ne "N0A0" ;
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
				
if($Info_hash{$filespecies}){
									
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




my $kmero;
my @trikmer1;
my $r = 'A,T,G,C';
for  my $kmer (glob"{$r}{$r}>{$r}{$r}"){
	## $kmer
	for ($kmer =~ m/(\w)(\w)>(\w)(\w)/xms){$kmero =join "",$1,'[',$2,'>',$3,']',$4};
	if ($2 ne $3){ push @trikmer1,$kmero };
}


my %total_mutation_context;
my %ID_mutation_count;
my $total_mut=0;




TAST:
for my $id (keys %result_hash){
	
	
	$ID_mutation_count{$id} = 0;
	#if( $id eq "N0A0" ){next TAST};
	my %res;
	my $time_data = $ID_times{$id};
	## $id
	my $seq = $result_hash{$id}{"seq"};
	my $nat = $result_hash{$id}{"Nature"};	
	my $io_ID = $id;
	$io_ID =~ s/N\d+//g;
	
	if ($id eq "N1A00"){
	
	$result_hash{$id}{"anc_seq"}=$result_hash{"N0A0"}{"seq"};
		## X:$result_hash{"N1A00"}
		
	};
	
	if ($id eq "N1A01"){
	
	$result_hash{$id}{"anc_seq"}=$result_hash{"N0A0"}{"seq"};
		## X:$result_hash{"N1A00"}
		
	};
	
	
	if ($result_hash{$id}{"anc_seq"}){
		
		## $id
		
		my $trueseq=$seq;
		$trueseq=~ s/\-//g;
		my $realseqlength = $the_seq_length; 
		for  my $trikmer1 (@trikmer1){
		
			$total_mutation_context{$io_ID}{$trikmer1}=0;
		
		};
		

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

	## $id

	@DNA_ref = split //, $seq_rf;

	my $statut = $result_hash{$id}{"Nature"};
	my $time_data = $time_values - $ID_times{$id};

	## $seq_rf
	## @DNA_ref





	TURN:
	while(@DNA_ref){
		my $nec = shift @DNA_ref;
		## $nec
		## $i
		
		if ($nec =~ m/\w/xms){ 
		for my $id (keys %{$res{$i}}){
			

			
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
					
					
					say {$out} join "\t",  $id,$nec, $i+2, $res{$i}{$id}{"Nucl"}, $res{$i}{$id}{"Codon"}, (join "",$up2,$up1,"[",$nec,">",$res{$i}{$id}{"Nucl"},"]",$down1,$down2),
					(join "",$nec,">",$res{$i}{$id}{"Nucl"}),$statut,$time_data,$full_time_values,$time_values,$realseqlength, $group,$genus,$family,$species,$filespecies;
					
					my $corrected_i=$i+1;
					
					$matrix_of_position{$corrected_i}++;
					
					my $mutation_type= join "", $up, "[", $nec, ">", $res{$i}{$id}{"Nucl"},"]",$down;
					
					
					
					$total_mutation_context{$io_ID}{$mutation_type}++;
					
					$ID_mutation_count{$id}++ ;
						
					
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




############################################################################################################################################################################################
				# %total_mutation_context
				
			## %ID_mutation_count



# %MATRIXES

my @values_kmers;


for  my $kmer (glob"{$r}C>T{$r}"){	
	
	my $CompI_kmero;	

	if ($kmer =~ m/(\w)(\w)\>(\w)(\w)/xms){
			
		$CompI_kmero =join "", $1,'[',$2,'>',$3,']',$4;
		
	};
	
	push @values_kmers,$CompI_kmero;

};

for  my $kmer (glob"{$r}C>G{$r}"){		
	
	my $CompI_kmero;	

	if ($kmer =~ m/(\w)(\w)\>(\w)(\w)/xms){
			
		$CompI_kmero =join "", $1,'[',$2,'>',$3,']',$4;
		
	};
	
	push @values_kmers,$CompI_kmero;

};

for  my $kmer (glob"{$r}C>A{$r}"){		
	
	my $CompI_kmero;	

	if ($kmer =~ m/(\w)(\w)\>(\w)(\w)/xms){
			
		$CompI_kmero =join "", $1,'[',$2,'>',$3,']',$4;
		
	};
	
	push @values_kmers,$CompI_kmero;

};


for  my $kmer (glob"{$r}T>A{$r}"){		
	
	my $CompI_kmero;	

	if ($kmer =~ m/(\w)(\w)\>(\w)(\w)/xms){
			
		$CompI_kmero =join "", $1,'[',$2,'>',$3,']',$4;
		
	};
	
	push @values_kmers,$CompI_kmero;

};


for  my $kmer (glob"{$r}T>G{$r}"){		
	
	my $CompI_kmero;	

	if ($kmer =~ m/(\w)(\w)\>(\w)(\w)/xms){
			
		$CompI_kmero =join "", $1,'[',$2,'>',$3,']',$4;
		
	};
	
	push @values_kmers,$CompI_kmero;

};


for  my $kmer (glob"{$r}A>G{$r}"){		
	
	my $CompI_kmero;	

	if ($kmer =~ m/(\w)(\w)\>(\w)(\w)/xms){
			
		$CompI_kmero =join "", $1,'[',$2,'>',$3,']',$4;
		
	};
	
	push @values_kmers,$CompI_kmero;

};

#################

my @listes_ifio= keys %total_mutation_context;


my %combined_muts;
my %Counting_combined_muts;

for  my $trikmer1 (@trikmer1){
	
	push @{ $combined_muts{$trikmer1}} , $trikmer1 ;
	
};


push @{ $combined_muts{"IDs"} }, "Mutation.Types" ;





for my $id (keys %total_mutation_context){
	
	my @seleced_IDs = grep (/${id}/, @listes_ifio);
	my %tmp_matrix;
		
	for  my $trikmer1 (@trikmer1){
	
		$tmp_matrix{$trikmer1}=0;
		$Counting_combined_muts{$id}{$trikmer1}=0;
	
	};
	
	$Counting_combined_muts{$id}{"TOTAL"}=0;
	
	for  my $seleced_IDs (@seleced_IDs){
		
		for  my $seleced_kmer ( keys %{$total_mutation_context{$seleced_IDs}} ){
		
			$tmp_matrix{$seleced_kmer}+= $total_mutation_context{$seleced_IDs}{$seleced_kmer};
			$Counting_combined_muts{$id}{$seleced_kmer}+= $total_mutation_context{$seleced_IDs}{$seleced_kmer};
			$Counting_combined_muts{$id}{"TOTAL"}+= $total_mutation_context{$seleced_IDs}{$seleced_kmer};
			
		};
		
	};
	
	push @{$combined_muts{"IDs"}}, $id ;
	
	for  my $trikmer1 (@trikmer1){
		
		push @{$combined_muts{$trikmer1}}, $tmp_matrix{$trikmer1};
		
	};
	

};

my $outname_final_mirror = join "", $outname, "SPLIT_Mirror_matrix.txt";
open my $out_final_mirror, '>',  $outname_final_mirror;

say {$out_final_mirror} join "\t",@{$combined_muts{"IDs"}};


for  my $trikmer1 (@trikmer1){

	say {$out_final_mirror} join "\t",@{$combined_muts{$trikmer1 }};

};

close $out_final_mirror;



######################################################################################################
my %combined_muts_prop;

my $outname_final_similar = join "", $outname, "_Score_similarity.txt";
open my $out_final_similar, '>',  $outname_final_similar;

say {$out_final_similar} join  "\t", "#OTU","Iteration","ID","True ID","Dist_ID","Cosine","Cosine_prop","Cosine_local","Nature","Branche_size","Branche_time","Tree_time","Measuring_timescale","branche_mut","branche_full_mut","total_mut",
									"Seq_length","Species","Group","Genus","Family";

my %sim;
my %similarity_prop;
my %similarity_simple;
##ROOT

my $id0="A0";
my $true_id0 = "N0A0";
my %tmp_matrix0;
my $tmp_tot0=0;
my @Dist0;





#my %petit_hash0 = $total_mutation_context{$id0};
	
#my @fw_values_simples0;
#my @rv_values_simples0;
	
#for  my $trikmer1 (@values_kmers){
		
#	my $Comp_kmero;
	
#	if ($trikmer1 =~ m/(\w)\[(\w)\>(\w)\](\w)/xms){
		
	#	$Comp_kmero =join "", $1,'[',$3,'>',$2,']',$4;
		
#	};
		
#	push @fw_values_simples0,$petit_hash0{$trikmer1};
#	push @rv_values_simples0,$petit_hash0{$Comp_kmero};
		
#};
	
#my $similarity_simple0 = get_cos_sim(96, \@fw_values_simples0, \@rv_values_simples0);

my $similarity_simple0 = 0;



		
for  my $trikmer1 (@trikmer1){
	
	$tmp_matrix0{$trikmer1}=0;
	
};
		
for  my $seleced_IDs (@listes_ifio){
		
	for  my $seleced_kmer ( keys %{$total_mutation_context{$seleced_IDs}} ){

		$tmp_tot0 += $total_mutation_context{$seleced_IDs}{$seleced_kmer};
		$tmp_matrix0{$seleced_kmer}+= $total_mutation_context{$seleced_IDs}{$seleced_kmer};
			
	};
		
};






	
my @fw_values0;
my @rv_values0;
	
for  my $trikmer1 (@values_kmers){
		
	my $Comp_kmero;
	
	if ($trikmer1 =~ m/(\w)\[(\w)\>(\w)\](\w)/xms){
			
		$Comp_kmero =join "", $1,'[',$3,'>',$2,']',$4;
		
	};
		
	my $prop = 0;
	my $prop_comp = 0;
		
	if ($tmp_tot0>0){
			
		$prop = ($tmp_matrix0{$trikmer1}/$tmp_tot0)*100;				
		$prop_comp = ($tmp_matrix0{$Comp_kmero}/$tmp_tot0)*100;
		
	};
		
	my $dist = abs ($prop - $prop_comp);
		
	push @Dist0,$dist ;
		
	push @fw_values0,$tmp_matrix0{$trikmer1};
	push @rv_values0,$tmp_matrix0{$Comp_kmero};
		
};
	
my $cosine0 = Bag::Similarity::Cosine->new;
my $similarity0 = get_cos_sim(96, \@fw_values0, \@rv_values0);
		
	############
my $DOWN10 = "A00";
my $DOWN20 = "A01";
		
my @fw_values_prop0;
my @rv_values_prop0;

my $tot_prop =0;
	
for  my $trikmer1 (@values_kmers){
		
	my $Comp_kmero;
	
	if ($trikmer1 =~ m/(\w)\[(\w)\>(\w)\](\w)/xms){
			
		$Comp_kmero =join "", $1,'[',$3,'>',$2,']',$4;
		
	};
		
	my $val_Down1 = $Counting_combined_muts{$DOWN10}{$trikmer1};

	my $val_Down2 = $Counting_combined_muts{$DOWN20}{$trikmer1};
	
	my $prop_Down1 = $val_Down1 /$Counting_combined_muts{$DOWN10}{"TOTAL"};
	my $prop_Down2 = $val_Down2 /$Counting_combined_muts{$DOWN20}{"TOTAL"};
	$tot_prop += $prop_Down1;
	$tot_prop += $prop_Down2;

	my $Final_val = ($prop_Down1 + $prop_Down2 )/2;
	
	push @fw_values_prop0, $Final_val ;
		
	my $compo_val_Down1 = $Counting_combined_muts{$DOWN10}{$Comp_kmero};
	my $compo_val_Down2 = $Counting_combined_muts{$DOWN20}{$Comp_kmero};

	my $compo_prop_Down1 = $compo_val_Down1 /$Counting_combined_muts{$DOWN10}{"TOTAL"};
	my $compo_prop_Down2 = $compo_val_Down2 /$Counting_combined_muts{$DOWN20}{"TOTAL"};
		
	my $compo_Final_val = ($compo_prop_Down1 + $compo_prop_Down2 )/2;
		
	push @rv_values_prop0,$compo_Final_val;
		
};


my $similarity_prop0 = get_cos_sim(96, \@fw_values_prop0, \@rv_values_prop0);
	 

	
	
#############
	
	
$similarity_simple{$true_id0}=$similarity_simple0;	
$similarity_prop{$true_id0}=$similarity_prop0;	
$sim{$true_id0}=$similarity0;

my $mpm0=sum(@Dist0);
my $final_dist0 = 100 - $mpm0 ;
		
say {$out_final_similar} join "\t",
$outname,
$iteratio,
$id0,
$true_id0,
$final_dist0,
$similarity0,
$similarity_prop0,
$similarity_simple0,	
$result_hash{$true_id0}{"Nature"},
$result_hash{$true_id0}{"Date"},
0,
$full_time_values,
$Measuring_timescale,
0,
$tmp_tot0,
$N_total,
$the_seq_length,
$species,
$group,
$genus,
$family;
		


TARE:
for my $id (keys %total_mutation_context){
	

	my $little_siez = (length $id)-2;
	my $true_id = join "", "N", $little_siez  ,$id;
	
	my %petit_hash = %{$total_mutation_context{$id}};
	
	
	my @fw_values_simples;
	my @rv_values_simples;
	
	for  my $trikmer1 (@values_kmers){
		
		my $Comp_kmero;
	
		if ($trikmer1 =~ m/(\w)\[(\w)\>(\w)\](\w)/xms){
			
			$Comp_kmero =join "", $1,'[',$3,'>',$2,']',$4;
		
		};
		
		push @fw_values_simples,$petit_hash{$trikmer1};
		push @rv_values_simples,$petit_hash{$Comp_kmero};
		
	};
	
	my $similarity_simple = get_cos_sim(96, \@fw_values_simples, \@rv_values_simples);
	
	if ($similarity_simple == 0 ){
		
		## $id
		## @fw_values_simples
		## @rv_values_simples
		
		## X:$result_hash{$true_id}{"anc_seq"}
		## y:$result_hash{$true_id}{"seq"}
	};
	
	##########
	
	
	my @seleced_IDs = grep (/${id}/, @listes_ifio);
	my $similarity_prop =0;
	my %tmp_matrix;
	my $tmp_tot=0;
	my @Dist;
		
	for  my $trikmer1 (@trikmer1){
	
		$tmp_matrix{$trikmer1}=0;
	
	};
		
	for  my $seleced_IDs (@seleced_IDs){
		
		for  my $seleced_kmer ( keys %{$total_mutation_context{$seleced_IDs}} ){
						## X:$tmp_matrix{$seleced_kmer}

			$tmp_tot += $total_mutation_context{$seleced_IDs}{$seleced_kmer};
			$tmp_matrix{$seleced_kmer}+= $total_mutation_context{$seleced_IDs}{$seleced_kmer};
			
		};
		
	};
	
	
	
	
	## %tmp_matrix
	
	my @fw_values;
	my @rv_values;
	
	for  my $trikmer1 (@values_kmers){
		
		my $Comp_kmero;
	
		if ($trikmer1 =~ m/(\w)\[(\w)\>(\w)\](\w)/xms){
			
			$Comp_kmero =join "", $1,'[',$3,'>',$2,']',$4;
		
		};
		
		my $prop = 0;
		my $prop_comp = 0;
		
		if ($tmp_tot>0){
			
			$prop = ($tmp_matrix{$trikmer1}/$tmp_tot)*100;				
			$prop_comp = ($tmp_matrix{$Comp_kmero}/$tmp_tot)*100;
		
		};
		
		my $dist = abs ($prop - $prop_comp);
		
		push @Dist,$dist ;
		
		push @fw_values,$tmp_matrix{$trikmer1};
		push @rv_values,$tmp_matrix{$Comp_kmero};
		
	};



	my $cosine = Bag::Similarity::Cosine->new;
	
	my $similarity = get_cos_sim(96, \@fw_values, \@rv_values);	

	
	############
	
	my $statut = $result_hash{$true_id}{"Nature"};
	## $statut 

	if ($statut =~ m/ncestor/xms){
		
		
		my $DOWN1 = join "",$id,0;
		my $DOWN2 = join "",$id,1;
			
		my @fw_values_prop;
		my @rv_values_prop;
		
		for  my $trikmer1 (@values_kmers){
			
			my $Comp_kmero;
		
			if ($trikmer1 =~ m/(\w)\[(\w)\>(\w)\](\w)/xms){
				
				$Comp_kmero =join "", $1,'[',$3,'>',$2,']',$4;
			
			};
			## $DOWN1
			## $trikmer1
			my $val_Down1 = $Counting_combined_muts{$DOWN1}{$trikmer1};
			my $val_Down2 = $Counting_combined_muts{$DOWN2}{$trikmer1};
			
			my $prop_Down1 ;

			if ( $val_Down1 == 0){
			
				$prop_Down1 = 0;
			
			}else{
				
				$prop_Down1 = $val_Down1 /$Counting_combined_muts{$DOWN1}{"TOTAL"};
			
			};
			
			my $prop_Down2 ;
			if ( $val_Down2 == 0){
				
				$prop_Down2 = 0;
				
			}else{
				
				$prop_Down2 = $val_Down2 /$Counting_combined_muts{$DOWN2}{"TOTAL"};
				
			};

			my $summ = $prop_Down1 + $prop_Down2;
			my $Final_val =0;
			if ( $summ > 0 ){
				
				$Final_val = ($prop_Down1 + $prop_Down2 )/2;
			
			};
			
			push @fw_values_prop, $Final_val ;
			
			my $compo_val_Down1 = $Counting_combined_muts{$DOWN1}{$Comp_kmero};
			my $compo_val_Down2 = $Counting_combined_muts{$DOWN2}{$Comp_kmero};
			
			my $compo_prop_Down1;
			my $compo_prop_Down2;
			
			if ( $compo_val_Down1 == 0){
			
				$compo_prop_Down1 = 0;
			
			}else{
				
				$compo_prop_Down1 = $compo_val_Down1 /$Counting_combined_muts{$DOWN1}{"TOTAL"};
			
			};
			

			if ( $compo_val_Down2 == 0){
				
				$compo_prop_Down2 = 0;
				
			}else{
				
				$compo_prop_Down2 = $compo_val_Down2 /$Counting_combined_muts{$DOWN2}{"TOTAL"};
				
			};

			my $compo_summ = $compo_prop_Down1 + $compo_prop_Down2;
			my $compo_Final_val =0;
			if ( $compo_summ > 0 ){
				
				$compo_Final_val = ($compo_prop_Down1 + $compo_prop_Down2 )/2;
			
			};
			
			push @rv_values_prop,$compo_Final_val;
			
		};

		my $test1 = sum @fw_values_prop;
		my $test2 = sum @rv_values_prop;		
		
		if ($test1 >0 and $test2>0){
		
			$similarity_prop = get_cos_sim(96, \@fw_values_prop, \@rv_values_prop);
		
		}else{
		
			$similarity_prop = 0;
		
		};
		
		$sim{$true_id}=$similarity;
		$similarity_prop{$true_id}=$similarity_prop;
		
	

	
	}else{
		
		$sim{$true_id}=0;
		$similarity_prop{$true_id}=0;
	
	};


	$similarity_simple{$true_id}=$similarity_simple;	
	
	if ($ID_mutation_count{$true_id}){
	
	## @Dist
		my $mpm=sum(@Dist);
		
		
		## $mpm
		my $final_dist = 100 - $mpm ;
		
		## $final_dist
		say {$out_final_similar} join "\t",
		$outname,
		$iteratio,
		$id,
		$true_id,
		$final_dist,
		$similarity,	
		$similarity_prop ,
		$similarity_simple,
		$result_hash{$true_id}{"Nature"},
		$result_hash{$true_id}{"Date"},
		$ID_times{$true_id},
		$full_time_values,
		$Measuring_timescale,
		$ID_mutation_count{$true_id},
		$tmp_tot,
		$N_total,
		$the_seq_length,
		$species,
		$group,
		$genus,
		$family;
		
	};

};



close $out_final_similar;











###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################

my $outname_final_similar2 = join "", $outname, "_Score_similarity2.txt";
open my $out_final_similar2, '>',  $outname_final_similar2;


say {$out_final_similar2} join  "\t", "#OTU","Iteration","ID","Dist_ID","Cosine","Nature","Branche_size","Branche_time","Tree_time","Measuring_timescale","branche_mut","branche_full_mut","total_mut",
									"Seq_length","Species","Group","Genus","Family";


		
say {$out_final_similar2} join "\t",
$outname,
$iteratio,
$id0,
$true_id0,
$similarity0,	
$result_hash{$true_id0}{"Nature"},
$result_hash{$true_id0}{"Date"},
0,
$full_time_values,
$Measuring_timescale,
0,
$tmp_tot0,
$N_total,
$the_seq_length,
$species,
$group,
$genus,
$family;
		


close $out_final_similar2;










##########################################################################################################################################################################
































###############################################



#my %total_mutation_context_mATRIX;
#my %tOTLA_ID_mutation_count;
#my @Secondelistes_ifio = keys %result_hash;

#TAST:
#for my $id (keys %result_hash){
	
#	my $statut = $result_hash{$id}{"Nature"};
	
	#if ($statut =~ m/ncestor/xms){
		
#		$tOTLA_ID_mutation_count{$id} =0;
#		for  my $trikmer1 (@trikmer1){
			
#			$total_mutation_context_mATRIX{$id}{$trikmer1}=0;
			
#		};
		
		
#		my $seq_rf = $result_hash{$id}{"seq"};
#		my @DNA_ref = split //, $seq_rf;
#		my ( $id_search ) = $id =~ m/(A\d+)/g;
#		my @sub_IDs = grep (/${id_search}/, @Secondelistes_ifio);
		
#		for my $subId (@sub_IDs){
			
#			my %res;
#			my $seq = $result_hash{$subId}{"seq"};
#			my $trueseq=$seq;
#			$trueseq=~ s/\-//g;
#			my $len = length $seq;
#			$seq=~ s/N/-/g;
	
#			PRE:
#			for(my $i = 0 ; $i <= $len ; $i += 1){
				
#				my $last = 0;
#				my $codon = substr ($seq, $i, 3);
				
#				if ($codon =~  m/\-(\w\w)/xms){
					
#					$codon=join "", $last, $1;
				
#				};

#				if ($codon =~  m/^\w/xms){
					
#					if ($codon =~  m/\w(\w)\-/xms){
					
#						$last = $1;
#						next PRE; 
					
#					}
								
#					if ($codon =~  m/\w(\w)\w/xms){
#							$res{$i}{$id} = {
#								"Nucl" => $1,
#								"Codon" => $codon,
#							};
							
#					};
					

#				};

#				next PRE;
				
#			};


#			my $i = -1;

#			## %res
#			my $statut = $result_hash{$id}{"Nature"};
		
#			TURN23:
#			while(@DNA_ref){
				 
#				my $nec = shift @DNA_ref;
				
#				if ($nec =~ m/\w/xms){ 
					
#					for my $id23 (keys %{$res{$i}}){

#						if($res{$i}{$id23}{"Nucl"} ne $nec){
							
#							my $up;
#							my $down;
			
#							for($res{$i}{$id23}{"Codon"} =~ m/(\w)\w(\w)/xms){$up = $1; $down = $2};  
									
#							my $mutation_type= join "", $up, "[", $nec, ">", $res{$i}{$id23}{"Nucl"},"]",$down;
			
#							$total_mutation_context_mATRIX{$id}{$mutation_type}++;
#							$tOTLA_ID_mutation_count{$id}++ ;

#							next;
							
#						};
						
#					};
						
#				};
				
#				$i = $i + 1;

#				next TURN23;
				
#			};

#		};
		
#	#};
	
#};


#my %TOTsim;

#TAREFULL:
#for my $id (keys %total_mutation_context_mATRIX){
	## $id 
#	my %tmp_matrix;
	# = $total_mutation_context_mATRIX{$id};
#	my $tmp_tot=0;
#	my @Dist;
#	my @fw_values;
#	my @rv_values;
	
#	for  my $trikmer1 (@values_kmers){
		
#		my $Comp_kmero;
	
#		if ($trikmer1 =~ m/(\w)\[(\w)\>(\w)\](\w)/xms){
			
#			$Comp_kmero =join "", $1,'[',$3,'>',$2,']',$4;
		
#		};
		
#		my $prop = 0;
#		my $prop_comp = 0;
		
#		if ($tmp_tot>0){
			
#			$prop = ($total_mutation_context_mATRIX{$id}{$trikmer1}/$tOTLA_ID_mutation_count{$id})*100;				
#			$prop_comp = ($total_mutation_context_mATRIX{$id}{$Comp_kmero}/$tOTLA_ID_mutation_count{$id})*100;
		
#		};
		
#		my $dist = abs ($prop - $prop_comp);	
#		push @Dist,$dist ;
#		push @fw_values,$total_mutation_context_mATRIX{$id}{$trikmer1};
#		push @rv_values,$total_mutation_context_mATRIX{$id}{$Comp_kmero};
		
#	};
	
#	my $cosine = Bag::Similarity::Cosine->new;	
	
#	## @fw_values
	## @rv_values
#	my $similarity = get_cos_sim(96, \@fw_values, \@rv_values);

#	$TOTsim{$id}=$similarity;

#};


################


#my $outname_final_full_mirror = join "", $outname, "Full_Mirror_matrix.txt";
#open my $out_final_mirror_tot, '>',  $outname_final_full_mirror;

#my @tmp_keys = keys %total_mutation_context_mATRIX;

#say {$out_final_mirror_tot} join "\t","Mutation.Types",@tmp_keys;



	
#for  my $trikmer1 (@trikmer1){
		
#	my @tmp_value;
#	for  my $tmpID2 (@tmp_keys){

#		push @tmp_value,$total_mutation_context_mATRIX{$tmpID2}{$trikmer1};
		
#	};
		
#	say {$out_final_mirror_tot} join "\t",$trikmer1 ,@tmp_value;
	
#};



#close $out_final_mirror_tot;



##########################################################################################################################################################################

## %sim


open my $in, '<', $tree_file;		  
my $tree_ref2 = <$in>;
chomp $tree_ref2;


my $YY=0;

my $outfile_tree = join "" , $outname,'_Zera.tree';
open my $out_tree, '>',  $outfile_tree;




say {$out_tree} "#NEXUS";
say {$out_tree} "";
say {$out_tree} "Begin taxa;";
say {$out_tree} join "","        Dimensions ntax=",scalar keys %result_hash ,";";
say {$out_tree} "        Taxlabels";




TUST:
for my $id (keys %result_hash){

	my $nat = $result_hash{$id}{"Nature"};
	
	say {$out_tree} join "","                ",$nat;
	
};

say {$out_tree} "                ;";
say {$out_tree} "End;";
say {$out_tree} "";
say {$out_tree} "Begin trees;";
#say {$out_tree} "[keywords: discretized_branch_rates]";
say {$out_tree} "        Translate";


## %similarity_prop

TUST:
for my $id (keys %result_hash){
	$YY++;
	
	## $id 
	## $YY
	#if( $id eq "N0A0" ){next TAST};
	my %res;
	my $time_data = $ID_times{$id};

	my $seq = $result_hash{$id}{"IDENT"};
	my $nat = $result_hash{$id}{"Nature"};
	
	## $nat
	$seq =~ s/\[/\\\[/g;
	$seq =~ s/\]/\\\]/g;
	$seq =~ s/\&/\\\&/g;
	$seq =~ s/\:/\\\:/g;
	
	$seq =~ s{\\\\}{\\}g;
	
	my $zim0 = 0;
	if ($sim{$id}){
		 $zim0 = $sim{$id};
	};
	
	my $zim1 = 0;
	if ($similarity_prop{$id}){
		 $zim1 = $similarity_prop{$id};
	};
	
	my $zim2 = 0;
	if ($similarity_simple{$id}){
		 $zim2 = $similarity_simple{$id};
	};
	
#	my $zim = 0;
#	if ($TOTsim{$id}){
		
#		 $zim = $TOTsim{$id};
	
#	};
	
	## $seq
	$tree_ref2 =~ s/(\d+)?$seq/$id\:\[\&cosine1\=$zim0\,\&cosine2\=$zim1\,\&cosine3\=$zim2\]/g ;
	## tata
	say {$out_tree} join "","		",$id,"\ ",$nat,",";
	
};

say {$out_tree} "                ;";


$tree_ref2 =~ s/\]\;/\]0\;/g;

say {$out_tree} $tree_ref2;

say {$out_tree} "End;";












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

	if ($tree_ref=~ m/.+(\[\&.+\=\"(.+)\"\]);/xms){


							
		my $date=0;
		my $veriRate =0;
		my $veri_seq = $2;
		my $identifier0 =$1;
		my $ANCEST0= join "", "Ancestor",$pop;
		$tree_ref =~ s/(.+)\[\&.+\=\"(.+)\"\];/$1/g;

		##  $veri_seq 

		$tree_hash{$e0}{$id01}= {
			"tree" => $tree_ref,
			"Value"=> "0",
			"seq"=> $veri_seq,
			"Nature" => $ANCEST0,
			"Date"=> $date,
			"Rate"=> $veriRate
		
		};
							
		$result_hash{$id01}= {
			"tree" => $tree_ref,
			"Value"=> "0",
			"seq"=> $veri_seq,
			#"anc_seq" => $tree_hash{"1"}{""}{"seq"},
			"Nature" => $ANCEST0,
			"Date"=> $date,
			"Rate"=> $veriRate,
			"IDENT"=> $identifier0
		
		};

							
	};
	
	
	## $tree_ref	
	## $Ngen	

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

				if($ext =~ m/\((.+)\)(\[&[^\=]+\=\"(\w+)\"\]\:\[.rate\=([^\]]+)\])([^\,^\)]+)/xms){
				
					my $date=lc($5);
					my $identifier = $2;
					my $veri_seq = $3;
					my $veriTree =$1;
					my $veriRate =lc ($4);
					
					my $ext_test  = $ext;
					$ext_test =~ s/[ATCG][ATCG][ATCG][ATCG][ATCG]+//g;
					
					## $ext_test 
					
					my $veriTree_test = $veriTree;
					$veriTree_test =~ s/[ATCG][ATCG][ATCG][ATCG][ATCG]+//g;
					
					## $veriTree_test
					
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
						
						my $In_mother_seq = $tree_hash{$y}{$id}{"seq"};
						
						
						my $ixt= $ext;
															
						$ixt =~ s/\(.+\)(\[&[^\=]+\=\"\w+\"\]\:\[.rate\=[^\]]+\][^\,^\)]+)/$1/g;
						
						my $axt= $ext;
					 
						$axt =~ s/\(.+\)(\[&[^\=]+\=\")\w+(\"\]\:\[.rate\=[^\]]+\][^\,^\)]+)/$1$new_id$2/g;
						
						if ( $id eq $id01 ){ $In_mother_seq =$veri_seq };
						
						my $mask= $veri_seq ^ $In_mother_seq;
						
						my @n_mut = ($mask =~ /[^\0]/g);
						
						my $n_mut = @n_mut;
						
						my $ANCEST= join "", "Ancestor",$pop;
						$pop++;


						$tree_hash{$e}{$new_id}= {
							"tree" => $veriTree,
							"seq"=> $veri_seq,
							"Nature"=>$ANCEST,
							"Date"=> $date,
							"Rate"=> $veriRate,
							"Nmute"=>$n_mut
						};
						
						$result_hash{$new_id}= {
							"anc_seq" => $tree_hash{$y}{$id}{"seq"},
							"seq"=> $veri_seq,
							"Nature"=>$ANCEST,
							"Date"=> $date,
							"Rate"=> $veriRate,
							"Nmute"=>$n_mut,
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
				

				
				if($ext =~ m/(\d+)(\[&[^\=]+\=\"(\w+)\"\]\:\[.rate\=([^\]]+)\])([^\,^\)]+)/xms){
					

					
					my $date=lc($5);
					my $identifier = $2;
					my $veri_seq = $3;
					my $numerous = $1;
					my $veriRate =lc ($4);
					
					my $ext_test  = $ext;
					$ext_test =~ s/[ATCG][ATCG][ATCG][ATCG][ATCG]+//g;
					
					## $ext_test 
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

## %result_hash
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

sub get_cos_sim {
	
    my ($size, $a, $b) = @_;
    my $n = $size-1;

    # Calculate eucledian magnitude
    my $a_len = 0;
    my $b_len = 0;
    foreach(0..$n) {
        $a_len += (@$a[$_]**2);
        $b_len += (@$b[$_]**2);
    }
    my $eucl_magn = sqrt($a_len * $b_len);
	my $dot_prod = 0;
    # If 0, stop calculation
    if ($eucl_magn == 0) {
        $dot_prod = 0;
    }else{

    # Calculate dot product
   
		foreach(0..$n) {
			$dot_prod += (@$a[$_] * @$b[$_]);
		}	
	}
    # Return cosine similarity
    my $sum=0;
    if ($eucl_magn> 0){
		
		$sum= $dot_prod / $eucl_magn;
	
	};
    
    return $sum;
    
}
