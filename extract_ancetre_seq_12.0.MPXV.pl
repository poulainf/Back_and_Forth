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
use Bag::Similarity::Cosine;

my $tree_file = shift;
my $ID_file = shift;

my $INFO_table= shift;
my $INFO_size= shift;

my $life_seq = shift;

##########

my $exclusion_file = shift;
open my $in_exclusion, '<', $exclusion_file;

my @exclusion;

LINE:
		
while (my $line_exclusion = <$in_exclusion> ){

	chomp $line_exclusion;

	push @exclusion, $line_exclusion;

};


#########
### mpx12

my ($outname) = $tree_file =~ m/(.*)\.\w+/g;
my ($outname_short) = "MPXV";

my ($filespecies) = $tree_file =~ m/(.*)\.fa.*\.tree*/g; 


my %Ancestor_hash = read_Anc( $ID_file );
read_Anc_life( $life_seq );


## x:keys(%Ancestor_hash)

## %ID_ref
my %Info_hash = read_info( $INFO_table);
my ($tree_hash, $result_hash, $ID_feuilles) = read_tree( $tree_file);
my %tree_hash = %$tree_hash;
my %result_hash = %$result_hash;
my %ID_feuilles = %$ID_feuilles; 


my %Info_length = read_length( $INFO_size);
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
my $N_total = 0;
my $full_time_values=1;
my $time_values = 1;
my $kj=0;
my $Measuring_timescale = 1;
		
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

my $selected_count=0;


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


### @exclusion

TAST:
for my $id (keys %result_hash){
	my $seq = $result_hash{$id}{"seq"};
	my $nat = $result_hash{$id}{"Nature"};	
	if ( grep /$nat/ , @exclusion){
		### ######################## STOP #################################################
		next TAST;
	};
	
	$ID_mutation_count{$id} = 0;
	if( $id eq "N0A0" ){next TAST};
	my %res;


	

	

	
	
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
	my $time_data = "NA";

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

my %Counting_muts;


for  my $trikmer1 (@trikmer1){
	
	push @{ $combined_muts{$trikmer1}} , $trikmer1 ;
	
	$Counting_muts{$trikmer1}=0;
	
};


push @{ $combined_muts{"IDs"} }, "Mutation.Types" ;

$Counting_muts{"TOTAL"}=0;



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
			$Counting_muts{$seleced_kmer}+= $total_mutation_context{$seleced_IDs}{$seleced_kmer};
			$Counting_muts{"TOTAL"}+= $total_mutation_context{$seleced_IDs}{$seleced_kmer};
			
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

my @list_ids0 = keys %Counting_combined_muts;


my $id0="A0";
my $true_id0 = "N0A0";
my %tmp_matrix0;
my $tmp_tot0=$Counting_muts{"TOTAL"};
my @Dist0;





my $similarity_simple0 = 0;



		
for  my $trikmer1 (@trikmer1){
	
	$tmp_matrix0{$trikmer1}=0;
	
};
		
for  my $seleced_IDs (@listes_ifio){
		
	for  my $seleced_kmer ( keys %{$total_mutation_context{$seleced_IDs}} ){

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
	
	my $val_Down = $Counting_muts{$trikmer1};
	my $prop_Down = $val_Down /$tmp_tot0;
	push @fw_values_prop0, $prop_Down ;
		
	my $compo_val_Down = $Counting_muts{$Comp_kmero};
	my $compo_prop_Down = $compo_val_Down /$tmp_tot0;	
	push @rv_values_prop0,$compo_prop_Down;
		
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
"1",
$id0,
$true_id0,
$final_dist0,
$similarity0,
$similarity_prop0,
$similarity_simple0,	
"ROOT",
1,
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
		"1",
		$id,
		$true_id,
		$final_dist,
		$similarity,	
		$similarity_prop ,
		$similarity_simple,
		$result_hash{$true_id}{"Nature"},
		1,
		0,
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
"1",
$id0,
$true_id0,
$similarity0,	
"ROOT",
1,
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


	my $TOTAL_A_count = 0;
	my $TOTAL_T_count = 0;
	my $TOTAL_C_count = 0;
	my $TOTAL_G_count = 0;
	my $pop =0;
	
	
	if($tree_ref =~ m/(.+)\)([^\)]*)\:[^\,]+?,[^\,]+\;/xms){
		
		### juju
		
		my $tree_ref = $1;
		my $ident=$2;
		my $veri_seq="NA";
		
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
		
		
		IDS:
		for my $id (keys %{$tree_hash{$y}}){
			

			my $tree = $tree_hash{$y}{$id}{"tree"};
			
			
			my $tree1= $tree;
			my $PATTERN = $RE{balanced}{-parens=>'()'};
			my @matches = ();
			
			INPUT:
			while (length $tree) {
				$tree =~ s{\A [^\w\(]+ }{}x;  ## skip leading junk
				if ($tree =~ s{ ($PATTERN[^\:]+\:[^\,^\)]+) }{}x ) {
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
		
	}else{
		
		## $a
		## $b
		
		## @$a
		## @$b
	
	
	};
    
    return $sum;
    
}
