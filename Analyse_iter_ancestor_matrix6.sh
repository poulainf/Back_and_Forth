#! /bin/bash
TREE_file=$1;
iter=$2;


TREE_file_true="$( echo $TREE_file |  grep -Po "[^\/]+$" )";
TREE_fixe="$( echo $TREE_file |  sed -e s/"${TREE_file_true}"// | sed -e s/"..\/"//)";
#TREE_fixe="./";
echo "$subname"


subname="$( echo $TREE_file_true |  sed -e s/.trees// )";
subname_true="$( echo $TREE_file |  sed -e s/.trees// | sed -e s/"..\/"//)";
#subname_true=$subname;
info_FILE=$3;

ID_files="${subname}_IDS.txt";
noBURNinTREES="${subname}_no_BURN.trees";

echo "Start of $TREE_file_true SNPs extraction ..."
date

echo "${TREE_fixe}"


echo $TREE_file 
	
	echo "no BURN-IN ..."
	line_start="$(grep "STATE" -Pno $TREE_file | head -n1 | cut -d ":" -f1)"
	((line_start=line_start-1))
	sed -e "1,$line_start d" $TREE_file > $noBURNinTREES;
	
	
grep -va "tree STATE" $TREE_file | grep -Poa "\d+\ [^\|]+\|[^\,]+"  > $ID_files;

for i in $(seq 1 $iter);do 	
	
	STATE="$(shuf -n 1 $noBURNinTREES | grep -Po "STATE_\d+")";
	TREE_files="Uniq_true_tree_${subname}_iteration_${i}.tree";
	Exclusion_file="${subname}Exclusion.txt";
	grep -a "$STATE " $noBURNinTREES   > $TREE_files;
	
	
	
	extract_ancetre_seq_9.0.pl $TREE_files $ID_files $info_FILE ../../Full_seq_length_line.txt;
	extract_ancetre_seq_12.0.pl $TREE_files $ID_files $info_FILE ../../Full_seq_length_line.txt;

done




#### Trinucl

NUM="$( ls Uniq_true_tree_${subname}_iteration_*matrix_line.txt | wc -l )";
y=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_line.txt `; do

	((y=y+1))
	tmp="tmp_${y}_"$subname"_trinuc.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_line.txt" > "Combine_${subname}_matrix.txt"; 
paste Combine_${subname}_matrix.txt tmp_*${subname}_trinuc.txt > "Final_Combine_${subname}_matrix.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix.txt" > "Final_Combine_${subname}_matrix2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix.txt" 

rm tmp_*"$subname"_trinuc.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_line.txt
rm Final_Combine_"$subname"_matrix2.txt
rm Combine_"$subname"_matrix.txt




#### Trinucl96

NUM96="$( ls Uniq_true_tree_${subname}_iteration_*matrix96_line.txt | wc -l )";
y96=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix96_line.txt `; do

	((y96=y96+1))
	tmp="tmp_${y96}_"$subname"_trinuc96.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix96_line.txt" > "Combine_${subname}_matrix96.txt"; 
paste Combine_${subname}_matrix96.txt tmp_*${subname}_trinuc96.txt > "Final_Combine_${subname}_matrix96.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix96.txt" > "Final_Combine_${subname}_matrix296.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix296.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix96.txt" 



rm tmp_*"$subname"_trinuc96.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix96_line.txt
rm Final_Combine_"$subname"_matrix296.txt
rm Combine_"$subname"_matrix96.txt


#### Penta

NUM_penta="$( ls Uniq_true_tree_${subname}_iteration_*matrix_penta_line.txt | wc -l )";
yp=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_penta_line.txt `; do

	((yp=yp+1))
	tmp="tmp_${yp}_"$subname"_penta.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_penta_line.txt" > "Combine_${subname}_matrix_penta_line.txt"; 
paste Combine_${subname}_matrix_penta_line.txt tmp_*${subname}_penta.txt > "Final_Combine_${subname}_matrix_penta_line.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_penta_line.txt" > "Final_Combine_${subname}_matrix_penta_line2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_penta_line2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_penta_line.txt" 



rm tmp_*"$subname"_penta.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_penta_line.txt
rm Final_Combine_"$subname"_matrix_penta_line2.txt
rm Combine_"$subname"_matrix_penta_line.txt



#### Single


NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_single.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_single.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_single.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_single.txt" > "Combine_${subname}_matrix_single.txt"; 
paste Combine_${subname}_matrix_single.txt tmp_*${subname}_single.txt > "Final_Combine_${subname}_matrix_single.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_single.txt" > "Final_Combine_${subname}_matrix_single2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_single2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_single.txt" 


rm tmp_*"$subname"_single.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_single.txt
rm Final_Combine_"$subname"_matrix_single2.txt
rm Combine_"$subname"_matrix_single.txt



### Double_nuc


NUM_double="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_double_line.txt | wc -l )";
yd6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_matrix_double_line.txt `; do

	((yd6=yd6+1))
	tmp="tmp_${yd6}_"$subname"_double.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_double_line.txt" > "Combine_${subname}_matrix_double_line.txt"; 
paste Combine_${subname}_matrix_double_line.txt tmp_*${subname}_double.txt > "Final_Combine_${subname}_matrix_double_line.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_double_line.txt" > "Final_Combine_${subname}_matrix_double_line2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_double_line2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_double_line.txt" 


rm tmp_*"$subname"_double.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_double_line.txt
rm Final_Combine_"$subname"_matrix_double_line2.txt
rm Combine_"$subname"_matrix_double_line.txt



#Uniq_true_tree_Fin_Corrected_Generated_Uniq_true_tree_NO-iT_parvovirus_b19_segment_taxon_10798_splited_iter-N0_Clustalo_ALIGNED_sub_selected_iteration__Ran_100__sub_selected_iteration_1_matrix_txE.txt





#### Sub rate
#Uniq_true_tree_NO-iT_hepatitis_b_virus_segment_subtype_ayw3_splited_iter-N0_Clustalo_ALIGNED_sub_selected_iteration_100_matrix_txS.txt

NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_txS.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_txS.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"TxS.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_txS.txt" > "Combine_${subname}_matrix_txS.txt"; 
paste Combine_${subname}_matrix_txS.txt tmp_*${subname}TxS.txt > "Final_Combine_${subname}_matrix_txS.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_txS.txt" > "Final_Combine_${subname}_matrix_txS2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_txS2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_txS.txt" 


rm tmp_*"$subname"TxS.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_txS.txt
rm Final_Combine_"$subname"_matrix_txS2.txt
rm Combine_"$subname"_matrix_txS.txt







NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_txR.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_txR.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"txR.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_txR.txt" > "Combine_${subname}_matrix_txR.txt"; 
paste Combine_${subname}_matrix_txR.txt tmp_*${subname}txR.txt > "Final_Combine_${subname}_matrix_txR.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_txR.txt" > "Final_Combine_${subname}_matrix_txR2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_txR2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_txR.txt" 


rm tmp_*"$subname"txR.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_txR.txt
rm Final_Combine_"$subname"_matrix_txR2.txt
rm Combine_"$subname"_matrix_txR.txt



NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_txE.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_txE.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"txE.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_txE.txt" > "Combine_${subname}_matrix_txE.txt"; 
paste Combine_${subname}_matrix_txE.txt tmp_*${subname}txE.txt > "Final_Combine_${subname}_matrix_txE.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_txE.txt" > "Final_Combine_${subname}_matrix_txE2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_txE2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_txE.txt" 


rm tmp_*"$subname"txE.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_txE.txt
rm Final_Combine_"$subname"_matrix_txE2.txt
rm Combine_"$subname"_matrix_txE.txt



# OTU COMPO
NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_compo.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_compo.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_compo.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_compo.txt" > "Combine_${subname}_matrix_compo.txt"; 
paste Combine_${subname}_matrix_compo.txt tmp_*${subname}_compo.txt > "Final_Combine_${subname}_matrix_compo.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_compo.txt" > "Final_Combine_${subname}_matrix_compo2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_compo2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_compo.txt" 


rm tmp_*"$subname"_compo.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_compo.txt
rm Final_Combine_"$subname"_matrix_compo2.txt
rm Combine_"$subname"_matrix_compo.txt








#### Sub rate pinta

NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_txS_penta.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_txS_penta.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_penta_txS.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_txS_penta.txt" > "Combine_${subname}_matrix_txS_penta.txt"; 
paste Combine_${subname}_matrix_txS_penta.txt tmp_*${subname}_penta_txS.txt > "Final_Combine_${subname}_matrix_txS_penta.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_txS_penta.txt" > "Final_Combine_${subname}_matrix_txS2_penta.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_txS2_penta.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_txS_penta.txt" 


rm tmp_*"$subname"_penta_txS.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_txS_penta.txt
rm Final_Combine_"$subname"_matrix_txS2_penta.txt
rm Combine_"$subname"_matrix_txS_penta.txt






NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_txR_penta.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_txR_penta.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_penta_txR.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

#rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_txR_penta.txt" > "Combine_${subname}_matrix_txR_penta.txt"; 
paste Combine_${subname}_matrix_txR_penta.txt tmp_*${subname}_penta_txR.txt > "Final_Combine_${subname}_matrix_txR_penta.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_txR_penta.txt" > "Final_Combine_${subname}_matrix_txR2_penta.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_txR2_penta.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_txR_penta.txt" 


rm tmp_*"$subname"_penta_txR.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_txR_penta.txt
rm Final_Combine_"$subname"_matrix_txR2_penta.txt
rm Combine_"$subname"_matrix_txR_penta.txt




NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_txE_penta.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_txE_penta.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_penta_txE.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_txE_penta.txt" > "Combine_${subname}_matrix_txE_penta.txt"; 
paste Combine_${subname}_matrix_txE_penta.txt tmp_*${subname}_penta_txE.txt > "Final_Combine_${subname}_matrix_txE_penta.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_txE_penta.txt" > "Final_Combine_${subname}_matrix_txE2_penta.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_txE2_penta.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_txE_penta.txt" 


rm tmp_*"$subname"_penta_txE.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_txE_penta.txt
rm Final_Combine_"$subname"_matrix_txE2_penta.txt
rm Combine_"$subname"_matrix_txE_penta.txt



# OTU COMPO pinta
NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_compo_penta.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_compo_penta.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_compo_pinta.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_compo_penta.txt" > "Combine_${subname}_matrix_compo_penta.txt"; 
paste Combine_${subname}_matrix_compo_penta.txt tmp_*${subname}_compo_pinta.txt > "Final_Combine_${subname}_matrix_compo_penta.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_compo_penta.txt" > "Final_Combine_${subname}_matrix_compo_pinta2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_compo_pinta2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_compo_penta.txt" 


rm tmp_*"$subname"_compo_pinta.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_compo_penta.txt
rm Final_Combine_"$subname"_matrix_compo_pinta2.txt
rm Combine_"$subname"_matrix_compo_penta.txt











#### Sub rate mono

NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_txS_mono.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_txS_mono.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_mono_txS.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_txS_mono.txt" > "Combine_${subname}_matrix_txS_mono.txt"; 
paste Combine_${subname}_matrix_txS_mono.txt tmp_*${subname}_mono_txS.txt > "Final_Combine_${subname}_matrix_txS_mono.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_txS_mono.txt" > "Final_Combine_${subname}_matrix_txS2_mono.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_txS2_mono.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_txS_mono.txt" 


rm tmp_*"$subname"_mono_txS.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_txS_mono.txt
rm Final_Combine_"$subname"_matrix_txS2_mono.txt
rm Combine_"$subname"_matrix_txS_mono.txt






NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_txR_mono.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_txR_mono.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_mono_txR.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_txR_mono.txt" > "Combine_${subname}_matrix_txR_mono.txt"; 
paste Combine_${subname}_matrix_txR_mono.txt tmp_*${subname}_mono_txR.txt > "Final_Combine_${subname}_matrix_txR_mono.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_txR_mono.txt" > "Final_Combine_${subname}_matrix_txR2_mono.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_txR2_mono.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_txR_mono.txt" 


rm tmp_*"$subname"_mono_txR.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_txR_mono.txt
rm Final_Combine_"$subname"_matrix_txR2_mono.txt
rm Combine_"$subname"_matrix_txR_mono.txt




NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*matrix_txE_mono.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_txE_mono.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_mono_txE.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_txE_mono.txt" > "Combine_${subname}_matrix_txE_mono.txt"; 
paste Combine_${subname}_matrix_txE_mono.txt tmp_*${subname}_mono_txE.txt > "Final_Combine_${subname}_matrix_txE_mono.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_txE_mono.txt" > "Final_Combine_${subname}_matrix_txE2_mono.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_txE2_mono.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_txE_mono.txt" 


rm tmp_*"$subname"_mono_txE.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_txE_mono.txt
rm Final_Combine_"$subname"_matrix_txE2_mono.txt
rm Combine_"$subname"_matrix_txE_mono.txt








# OTU COMPO mono

NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_compo_mono.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*matrix_compo_mono.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_compo_mona.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_compo_mono.txt" > "Combine_${subname}_matrix_compo_mono.txt"; 
paste Combine_${subname}_matrix_compo_mono.txt tmp_*${subname}_compo_mona.txt > "Final_Combine_${subname}_matrix_compo_mono.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_compo_mono.txt" > "Final_Combine_${subname}_matrix_compo_mono2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_compo_mono2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_compo_mono.txt" 


rm tmp_*"$subname"_compo_mona.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_compo_mono.txt
rm Final_Combine_"$subname"_matrix_compo_mono2.txt
rm Combine_"$subname"_matrix_compo_mono.txt





### nopingpong


NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*_position_rate.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_position_rate.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"position_rate.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_position_rate.txt" > "Combine_${subname}_position_rate.txt"; 
paste Combine_${subname}_position_rate.txt tmp_*${subname}position_rate.txt > "Final_Combine_${subname}_position_rate.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_position_rate.txt" > "Final_Combine_${subname}_position_rate2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_position_rate2.txt" > "${TREE_fixe}Final_Combine_${subname}_position_rate.txt" 


rm tmp_*"$subname"position_rate.txt
rm Uniq_true_tree_"$subname"_iteration_*_position_rate.txt
rm Final_Combine_"$subname"_position_rate2.txt
rm Combine_"$subname"_position_rate.txt






### pingpong


NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_192_PING-PONG_full.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_matrix_192_PING-PONG_full.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"matrix_192_PING-PONG_full.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_192_PING-PONG_full.txt" > "Combine_${subname}_matrix_192_PING-PONG_full.txt"; 
paste Combine_${subname}_matrix_192_PING-PONG_full.txt tmp_*${subname}matrix_192_PING-PONG_full.txt > "Final_Combine_${subname}_matrix_192_PING-PONG_full.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_192_PING-PONG_full.txt" > "Final_Combine_${subname}_matrix_192_PING-PONG_full2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_192_PING-PONG_full2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_192_PING-PONG_full.txt" 


rm tmp_*"$subname"matrix_192_PING-PONG_full.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_192_PING-PONG_full.txt
rm Final_Combine_"$subname"_matrix_192_PING-PONG_full2.txt
rm Combine_"$subname"_matrix_192_PING-PONG_full.txt





### nopingpong


NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_192_NO_PING-PONG_full.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_matrix_192_NO_PING-PONG_full.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"matrix_192_NO_PING-PONG_full.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_192_NO_PING-PONG_full.txt" > "Combine_${subname}_matrix_192_NO_PING-PONG_full.txt"; 
paste Combine_${subname}_matrix_192_NO_PING-PONG_full.txt tmp_*${subname}matrix_192_NO_PING-PONG_full.txt > "Final_Combine_${subname}_matrix_192_NO_PING-PONG_full.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_192_NO_PING-PONG_full.txt" > "Final_Combine_${subname}_matrix_192_NO_PING-PONG_full2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_192_NO_PING-PONG_full2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_192_NO_PING-PONG_full.txt" 


rm tmp_*"$subname"matrix_192_NO_PING-PONG_full.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_192_NO_PING-PONG_full.txt
rm Final_Combine_"$subname"_matrix_192_NO_PING-PONG_full2.txt
rm Combine_"$subname"_matrix_192_NO_PING-PONG_full.txt











########################### Matrix for deconstruction ##################################################################################



#### Sub rate_corr
NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txS_proportion.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txS_proportion.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"TxS_CORR.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_Corrected_matrix_txS_proportion.txt" > "Combine_${subname}_Corrected_matrix_txS_proportion.txt"; 
paste Combine_${subname}_Corrected_matrix_txS_proportion.txt tmp_*${subname}TxS_CORR.txt > "Final_Combine_${subname}_Corrected_matrix_txS_proportion.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_Corrected_matrix_txS_proportion.txt" > "Final_Combine_${subname}_Corrected_matrix_txS_proportion2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_Corrected_matrix_txS_proportion2.txt" > "${TREE_fixe}Final_Combine_${subname}_Corrected_matrix_txS_proportion.txt" 


rm tmp_*"$subname"TxS_CORR.txt
rm Uniq_true_tree_"$subname"_iteration_*_Corrected_matrix_txS_proportion.txt
rm Final_Combine_"$subname"_Corrected_matrix_txS_proportion2.txt
rm Combine_"$subname"_Corrected_matrix_txS_proportion.txt







NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txR_proportion.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txR_proportion.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"txR_CORR.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_Corrected_matrix_txR_proportion.txt" > "Combine_${subname}_Corrected_matrix_txR_proportion.txt"; 
paste Combine_${subname}_Corrected_matrix_txR_proportion.txt tmp_*${subname}txR_CORR.txt > "Final_Combine_${subname}_Corrected_matrix_txR_proportion.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_Corrected_matrix_txR_proportion.txt" > "Final_Combine_${subname}_Corrected_matrix_txR_proportion2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_Corrected_matrix_txR_proportion2.txt" > "${TREE_fixe}Final_Combine_${subname}_Corrected_matrix_txR_proportion.txt" 


rm tmp_*"$subname"txR_CORR.txt
rm Uniq_true_tree_"$subname"_iteration_*_Corrected_matrix_txR_proportion.txt
rm Final_Combine_"$subname"_Corrected_matrix_txR_proportion2.txt
rm Combine_"$subname"_Corrected_matrix_txR_proportion.txt



NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txE_proportion.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txE_proportion.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"txE_CORR.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_Corrected_matrix_txE_proportion.txt" > "Combine_${subname}_Corrected_matrix_txE_proportion.txt"; 
paste Combine_${subname}_Corrected_matrix_txE_proportion.txt tmp_*${subname}txE_CORR.txt > "Final_Combine_${subname}_Corrected_matrix_txE_proportion.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_Corrected_matrix_txE_proportion.txt" > "Final_Combine_${subname}_Corrected_matrix_txE_proportion2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_Corrected_matrix_txE_proportion2.txt" > "${TREE_fixe}Final_Combine_${subname}_Corrected_matrix_txE_proportion.txt" 


rm tmp_*"$subname"txE_CORR.txt
rm Uniq_true_tree_"$subname"_iteration_*_Corrected_matrix_txE_proportion.txt
rm Final_Combine_"$subname"_Corrected_matrix_txE_proportion2.txt
rm Combine_"$subname"_Corrected_matrix_txE_proportion.txt











#### Sub rate pinta

NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txS_penta_proportion.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txS_penta_proportion.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_penta_CORR_txS.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_Corrected_matrix_txS_penta_proportion.txt" > "Combine_${subname}_Corrected_matrix_txS_penta_proportion.txt"; 
paste Combine_${subname}_Corrected_matrix_txS_penta_proportion.txt tmp_*${subname}_penta_CORR_txS.txt > "Final_Combine_${subname}_Corrected_matrix_txS_penta_proportion.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_Corrected_matrix_txS_penta_proportion.txt" > "Final_Combine_${subname}_Corrected_matrix_txS2_penta_proportion.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_Corrected_matrix_txS2_penta_proportion.txt" > "${TREE_fixe}Final_Combine_${subname}_Corrected_matrix_txS_penta_proportion.txt" 


rm tmp_*"$subname"_penta_CORR_txS.txt
rm Uniq_true_tree_"$subname"_iteration_*_Corrected_matrix_txS_penta_proportion.txt
rm Final_Combine_"$subname"_Corrected_matrix_txS2_penta_proportion.txt
rm Combine_"$subname"_Corrected_matrix_txS_penta_proportion.txt






NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txR_penta_proportion.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txR_penta_proportion.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_penta_CORR_txR.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_Corrected_matrix_txR_penta_proportion.txt" > "Combine_${subname}_Corrected_matrix_txR_penta_proportion.txt"; 
paste Combine_${subname}_Corrected_matrix_txR_penta_proportion.txt tmp_*${subname}_penta_CORR_txR.txt > "Final_Combine_${subname}_Corrected_matrix_txR_penta_proportion.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_Corrected_matrix_txR_penta_proportion.txt" > "Final_Combine_${subname}_Corrected_matrix_txR2_penta_proportion.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_Corrected_matrix_txR2_penta_proportion.txt" > "${TREE_fixe}Final_Combine_${subname}_Corrected_matrix_txR_penta_proportion.txt" 


rm tmp_*"$subname"_penta_CORR_txR.txt
rm Uniq_true_tree_"$subname"_iteration_*_Corrected_matrix_txR_penta_proportion.txt
rm Final_Combine_"$subname"_Corrected_matrix_txR2_penta_proportion.txt
rm Combine_"$subname"_Corrected_matrix_txR_penta_proportion.txt




NUM_single="$( ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txE_penta_proportion.txt | wc -l )";
ys6=1

for i in `ls Uniq_true_tree_${subname}_iteration_*Corrected_matrix_txE_penta_proportion.txt `; do

	((ys6=ys6+1))
	tmp="tmp_${ys6}_"$subname"_penta_CORR_txE.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_Corrected_matrix_txE_penta_proportion.txt" > "Combine_${subname}_Corrected_matrix_txE_penta_proportion.txt"; 
paste Combine_${subname}_Corrected_matrix_txE_penta_proportion.txt tmp_*${subname}_penta_CORR_txE.txt > "Final_Combine_${subname}_Corrected_matrix_txE_penta_proportion.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_Corrected_matrix_txE_penta_proportion.txt" > "Final_Combine_${subname}_Corrected_matrix_txE2_penta_proportion.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_Corrected_matrix_txE2_penta_proportion.txt" > "${TREE_fixe}Final_Combine_${subname}_Corrected_matrix_txE_penta_proportion.txt" 


rm tmp_*"$subname"_penta_CORR_txE.txt
rm Uniq_true_tree_"$subname"_iteration_*_Corrected_matrix_txE_penta_proportion.txt
rm Final_Combine_"$subname"_Corrected_matrix_txE2_penta_proportion.txt
rm Combine_"$subname"_Corrected_matrix_txE_penta_proportion.txt









############## mATRIX NORMALIZED ################################

#### Trinucl

NUM="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_normalized_line.txt | wc -l )";
y=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_matrix_normalized_line.txt `; do

	((y=y+1))
	tmp="tmp_${y}_"$subname"_trinuc_norm.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_normalized_line.txt" > "Combine_${subname}_matrix_normalized.txt"; 
paste Combine_${subname}_matrix_normalized.txt tmp_*${subname}_trinuc_norm.txt > "Final_Combine_${subname}_matrix_normalized.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_normalized.txt" > "Final_Combine_${subname}_matrix_normalized2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_normalized2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_normalized.txt" 

rm tmp_*"$subname"_trinuc_norm.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_normalized_line.txt
rm Final_Combine_"$subname"_matrix_normalized2.txt
rm Combine_"$subname"_matrix_normalized.txt





NUM="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_normalized_positiv_line.txt | wc -l )";
y=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_matrix_normalized_positiv_line.txt `; do

	((y=y+1))
	tmp="tmp_${y}_"$subname"_trinuc_norm_plus.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_normalized_positiv_line.txt" > "Combine_${subname}_matrix_normalized_positiv.txt"; 
paste Combine_${subname}_matrix_normalized_positiv.txt tmp_*${subname}_trinuc_norm_plus.txt > "Final_Combine_${subname}_matrix_normalized_positiv.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_normalized_positiv.txt" > "Final_Combine_${subname}_matrix_normalized_positiv2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_normalized_positiv2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_normalized_positiv.txt" 

rm tmp_*"$subname"_trinuc_norm_plus.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_normalized_positiv_line.txt
rm Final_Combine_"$subname"_matrix_normalized_positiv2.txt
rm Combine_"$subname"_matrix_normalized_positiv.txt








NUM="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_normalized_negativ_line.txt | wc -l )";
y=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_matrix_normalized_negativ_line.txt `; do

	((y=y+1))
	tmp="tmp_${y}_"$subname"_trinuc_norm_neg.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_normalized_negativ_line.txt" > "Combine_${subname}_matrix_normalized_negativ.txt"; 
paste Combine_${subname}_matrix_normalized_negativ.txt tmp_*${subname}_trinuc_norm_neg.txt > "Final_Combine_${subname}_matrix_normalized_negativ.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_normalized_negativ.txt" > "Final_Combine_${subname}_matrix_normalized_negativ2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_normalized_negativ2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_normalized_negativ.txt" 

rm tmp_*"$subname"_trinuc_norm_neg.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_normalized_negativ_line.txt
rm Final_Combine_"$subname"_matrix_normalized_negativ2.txt
rm Combine_"$subname"_matrix_normalized_negativ.txt




#### Penta

NUM_penta="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_pinta_normalized_3086_line.txt | wc -l )";
yp=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_matrix_pinta_normalized_3086_line.txt `; do

	((yp=yp+1))
	tmp="tmp_${yp}_"$subname"_penta_norm.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_pinta_normalized_3086_line.txt" > "Combine_${subname}_matrix_pinta_normalized_3086_line.txt"; 
paste Combine_${subname}_matrix_pinta_normalized_3086_line.txt tmp_*${subname}_penta_norm.txt > "Final_Combine_${subname}_matrix_pinta_normalized_3086_line.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_pinta_normalized_3086_line.txt" > "Final_Combine_${subname}_matrix_pinta_normalized_3086_line2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_pinta_normalized_3086_line2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_pinta_normalized_3086_line.txt" 



rm tmp_*"$subname"_penta_norm.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_pinta_normalized_3086_line.txt
rm Final_Combine_"$subname"_matrix_pinta_normalized_3086_line2.txt
rm Combine_"$subname"_matrix_pinta_normalized_3086_line.txt




NUM_penta="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_pinta_normalized_3086_positiv_line.txt | wc -l )";
yp=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_matrix_pinta_normalized_3086_positiv_line.txt `; do

	((yp=yp+1))
	tmp="tmp_${yp}_"$subname"_penta_norm_positiv.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_pinta_normalized_3086_positiv_line.txt" > "Combine_${subname}_matrix_pinta_normalized_3086_positiv_line.txt"; 
paste Combine_${subname}_matrix_pinta_normalized_3086_positiv_line.txt tmp_*${subname}_penta_norm_positiv.txt > "Final_Combine_${subname}_matrix_pinta_normalized_3086_positiv_line.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_pinta_normalized_3086_positiv_line.txt" > "Final_Combine_${subname}_matrix_pinta_normalized_3086_positiv_line2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_pinta_normalized_3086_positiv_line2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_pinta_normalized_3086_positiv_line.txt" 



rm tmp_*"$subname"_penta_norm_positiv.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_pinta_normalized_3086_positiv_line.txt
rm Final_Combine_"$subname"_matrix_pinta_normalized_3086_positiv_line2.txt
rm Combine_"$subname"_matrix_pinta_normalized_3086_positiv_line.txt




NUM_penta="$( ls Uniq_true_tree_${subname}_iteration_*_matrix_pinta_normalized_3086_negativ_line.txt | wc -l )";
yp=1

for i in `ls Uniq_true_tree_${subname}_iteration_*_matrix_pinta_normalized_3086_negativ_line.txt `; do

	((yp=yp+1))
	tmp="tmp_${yp}_"$subname"_penta_norm_negativ.txt";
	rr="${tmp}p";
	cut -f2 $i > $tmp ; 

done

##rm tmp_1_"$subname".txt
cut -f1 "Uniq_true_tree_${subname}_iteration_1_matrix_pinta_normalized_3086_negativ_line.txt" > "Combine_${subname}_matrix_pinta_normalized_3086_negativ_line.txt"; 
paste Combine_${subname}_matrix_pinta_normalized_3086_negativ_line.txt tmp_*${subname}_penta_norm_negativ.txt > "Final_Combine_${subname}_matrix_pinta_normalized_3086_negativ_line.txt"; 
sed -e s/"_Clustalo_ALIGNED"//g "Final_Combine_${subname}_matrix_pinta_normalized_3086_negativ_line.txt" > "Final_Combine_${subname}_matrix_pinta_normalized_3086_negativ_line2.txt" 
sed -e s/"Uniq_true_tree_"//g "Final_Combine_${subname}_matrix_pinta_normalized_3086_negativ_line2.txt" > "${TREE_fixe}Final_Combine_${subname}_matrix_pinta_normalized_3086_negativ_line.txt" 



rm tmp_*"$subname"_penta_norm_negativ.txt
rm Uniq_true_tree_"$subname"_iteration_*_matrix_pinta_normalized_3086_negativ_line.txt
rm Final_Combine_"$subname"_matrix_pinta_normalized_3086_negativ_line2.txt
rm Combine_"$subname"_matrix_pinta_normalized_3086_negativ_line.txt











### Last_formatage

report_CONSENSUS="${subname_true}_mut-report.txt"
Title="$(head -n1 Uniq_true_tree_${subname}_iteration_1_mut-report_CONSENSUS.txt)"
echo "#File	$Title"  > $report_CONSENSUS
for lolo in ` ls Uniq_true_tree_${subname}_iteration_*_mut-report_CONSENSUS.txt `; do  sed -i 1d  $lolo ;awk -F '\t', '{print FILENAME ,"\t", $0}' $lolo  >>  $report_CONSENSUS ;done
rm Uniq_true_tree_"$subname"_iteration_*_mut-report_CONSENSUS.txt


report_fasta="${subname_true}_mut.fasta"
echo -n > $report_fasta
for y in `ls Uniq_true_tree_${subname}_iteration_*_seq.fasta`; do Iter="$( echo $y | grep -Po "iteration\_\d+" )"; cat $y | sed -e s/">"/">"$Iter"_"/ >> $report_fasta ; done
rm Uniq_true_tree_"$subname"_iteration_*_seq.fasta

report_tree="${subname_true}_mut.tree"
echo -n > $report_tree
for lolo in ` ls Uniq_true_tree_${subname}_iteration_*.tree`; do awk -F '\t', '{print FILENAME , $0}' $lolo >>  $report_tree ;done
# Uniq_true_tree_"$subname"_iteration_*.tree


report_DOUBLE="${subname_true}_mut-report_DOUBLE.txt"
Title="$( head -n1 Uniq_true_tree_${subname}_iteration_1_mut-report_DOUBLE.txt)"
echo "#File	$Title" > $report_DOUBLE
for lolo in ` ls Uniq_true_tree_${subname}_iteration_*_mut-report_DOUBLE.txt `; do  sed -i 1d  $lolo;  awk -F '\t', '{print FILENAME ,"\t", $0}' $lolo  >>  $report_DOUBLE ;done
rm Uniq_true_tree_"$subname"_iteration_*_mut-report_DOUBLE.txt




report_RATE="${subname_true}_Comparaison_sub_rate.txt"
Title="$( head -n1 Uniq_true_tree_${subname}_iteration_1_Comparaison_sub_rate.txt)"
echo "#File	$Title" > $report_RATE
for lolo in ` ls Uniq_true_tree_${subname}_iteration_*_Comparaison_sub_rate.txt `; do  sed -i 1d  $lolo;  awk -F '\t', '{print FILENAME ,"\t", $0}' $lolo  >>  $report_RATE ;done
rm Uniq_true_tree_"$subname"_iteration_*_Comparaison_sub_rate.txt


report_mono="${subname_true}_compilation_single.txt"
Title="$( head -n1 Uniq_true_tree_${subname}_iteration_1_compilation_single.txt)"
echo "#File	$Title" > $report_mono
for lolo in ` ls Uniq_true_tree_${subname}_iteration_*_compilation_single.txt `; do  sed -i 1d  $lolo;  awk -F '\t', '{print FILENAME ,"\t", $0}' $lolo  >>  $report_mono ;done
rm Uniq_true_tree_"$subname"_iteration_*_compilation_single.txt


report_mono_TxS="${subname_true}_compilation_TxS_mono.txt"
Title="$( head -n1 Uniq_true_tree_${subname}_iteration_1_compilation_TxS_mono.txt)"
echo "#File	$Title" > $report_mono_TxS
for lolo in ` ls Uniq_true_tree_${subname}_iteration_*_compilation_TxS_mono.txt `; do  sed -i 1d  $lolo;  awk -F '\t', '{print FILENAME ,"\t", $0}' $lolo  >>  $report_mono_TxS ;done
rm Uniq_true_tree_"$subname"_iteration_*_compilation_TxS_mono.txt



report_compile_sim2="${subname_true}_compilation_sim_scores2.txt"
head -n1 Uniq_true_tree_${subname}_iteration_1_Score_similarity2.txt > $report_compile_sim2
for lolo in ` ls Uniq_true_tree_${subname}_iteration_*_Score_similarity2.txt `; do  tail -n+2 $lolo >>  $report_compile_sim2 ;done
#rm Uniq_true_tree_"$subname"_iteration_*_Score_similarity2.txt


report_POSITION_PING="${subname_true}_matrix_POSITION_PING-PONG_full.txt"
head -n1 Uniq_true_tree_${subname}_iteration_1_matrix_POSITION_PING-PONG_full.txt > $report_POSITION_PING
for lolo in ` ls Uniq_true_tree_${subname}_iteration_*_matrix_POSITION_PING-PONG_full.txt `; do  tail -n+2 $lolo >>  $report_POSITION_PING ;done
rm Uniq_true_tree_"$subname"_iteration_*_matrix_POSITION_PING-PONG_full.txt

report_compile_sim="${subname_true}_compilation_sim_scores.txt"
head -n1 Uniq_true_tree_${subname}_iteration_1_Score_similarity.txt > $report_compile_sim
for lolo in ` ls Uniq_true_tree_${subname}_iteration_*_Score_similarity.txt `; do  tail -n+2 $lolo >>  $report_compile_sim ;done
#rm Uniq_true_tree_"$subname"_iteration_*_Score_similarity.txt


rm $noBURNinTREES
rm $ID_files

echo "Finished"
