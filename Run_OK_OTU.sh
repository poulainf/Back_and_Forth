#!/bin/bash

n=0
x=0
file="$1"
subname="$(echo $file | sed -e 's/.fa//')"

if [ -n "$2" ]; then 
    t="$2"
else 
    t="Cl"
fi

while [ $x -le 2 ]; do 
    n=$((n + 1))
    
    echo "Start splitting $subname iteration: N$n"
    Splite_Iter_fasta_files.pl "$file" 50 "$n"
    
    if [ "$t" = "Ka" ]; then
        echo "Start MSA $subname iteration: N$n / Method: Kalign"
        kalign "${subname}_splited_iter-N${n}.fa" -f fasta > "${subname}_splited_iter-N${n}tmp.fast"
        tail -n +27 "${subname}_splited_iter-N${n}tmp.fast" > "${subname}_splited_iter-N${n}_ALIGNED.fasta"
    else
        echo "Start MSA $subname iteration: N$n / Method: Clustalo"
        clustalo -i "${subname}_splited_iter-N${n}.fa" -o "${subname}_splited_iter-N${n}_ALIGNED.fasta"
    fi
    
    echo "Start Tree reconstruction $subname iteration: N$n"
    FastTree < "${subname}_splited_iter-N${n}_ALIGNED.fasta" > "${subname}_splited_iter-N${n}_ALIGNED.tree"
    
    echo "Start rerooting $subname iteration: N$n"
    ~/Genomics/stat/reroot.pl -midpoint < "${subname}_splited_iter-N${n}_ALIGNED.tree" > "Re_root_${subname}_splited_iter-N${n}_ALIGNED.tree"
    
    echo "Start p-value testing $subname iteration: N$n"
    /home/florian1/Desktop/analyze/Trie_fichier-beastOK.R "Re_root_${subname}_splited_iter-N${n}_ALIGNED.tree"
    
    if [ -f "Ok_${subname}_splited_iter-N${n}_ALIGNED_tree.txt" ]; then
        x=$((x + 1))
    fi
done
