# Virus Substitutions Collection Pipeline

## Generate Neighbor-to-RefSeq Sequence File 

In our investigation, we only selected sequences for human host viruses. The following steps outline how to generate the required file:

```bash
# Download the genome data file from NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Viruses_RefSeq_and_neighbors_genome_data.tab

# Filter the data to include only sequences associated with human hosts
grep human Viruses_RefSeq_and_neighbors_genome_data.tab > Human_Viruses_RefSeq_and_neighbors_genome_data.tab

# Clean and process the filtered data
Clean_Neighbor_Refseq_table.pl Human_Viruses_RefSeq_and_neighbors_genome_data.tab
```
These commands will produce a file named `Human_Viruses_RefSeq_and_neighbors_genome_data-cleaned.tsv`. This file associates each full genome sequence of a human virus from the NCBI nucleotide database with a specific reference sequence.

## OTU fasta files generation
The cleaned table file can then be used to download FASTA and GenBank files using `esearch` from the Entrez [Direct utilitiez](https://www.ncbi.nlm.nih.gov/books/NBK179288/).The following command allow an automated files download and concatenate them. Sequence genomes clustering into Operative Taxonomic Units are next performed based on taxonomic information provided by GeneBank datas and Neighbor-to-RefSeq Sequence File. 

```bash
# Download the selected genome data files and GeneBank files from NCBI. The package produce one directory per OTU.
Download_Fasta_GeneBank.pl Human_Viruses_RefSeq_and_neighbors_genome_data-cleaned.tsv

# Concatenate files
for dir in */; do
    dir_name="${dir%/}"
    cat "$dir_name"/*.fasta > "${dir_name}.fasta"
done

# Cluster sequence into OTU
Ref_Fasta_to_OTU.pl test.fasta test.gb test.tsv
```
This step result on the production of large number of OTUs fasta files contening all the fasta sequences corresponding to a virus subspecies taxonomical scale. 

## OTU Selection for Temporal Signal Testing

All OTU full FASTA files are used to generate sub-FASTA files containing a maximum of 50 genome sequences, which are tested for temporal signal and prepared for phylogenetic inference. For OTU full FASTA files with more than 50 sequences, three randomly selected sub-FASTA files are generated per OTU.
Requirements

    Software Installation:
    Ensure the local installation of the following tools, with their executables accessible via the environment path:
        Clustal Omega
        Kalign
        FastTree

    Script Dependencies:
    The execution of the OTU_split_and_temporal_signal_test.sh script requires in the current directory:
        The Genomics/ and Temporal_signal_functions/ directories.
        The Trie_fichier-beastOK.R script.

```bash
# Run split OTU
for i in `ls *fasta`;do

  OTU_split_and_temporal_signal_test.sh $i ;

done
```

## OTUs phylogenetic trees inference
To perform BEAST1 phylogenetic tree inferences automatically, XML files are generated. Additionally, a Bash script is created for each OTU to enable the execution of the package in a Slurm-based cluster environment.

```bash
for i in ` ls *_ALIGNED.fasta `; do

  beauty_date_constantpop_GTR+G+I.pl $i 1000000000000 100000;

done

for i in `ls *.xml | sed -e s/"_formated.xml"//`; do

  ((a = a + 1));
  sub="$( echo $i | cut -d "_" -f2 | cut -c1-5 )";
  id="${sub}_${a}" ;
  tpage --define name=$i --define nick=$id creat_run.tt > Submition_$i.sh;

done
```
Burn-in values are determined through manual inspection using the Tracer package and recorded in the Burning.txt file.
Based on the Burning.txt file, the file.trees is trimmed for burn-in, and 500 trees are selected.

```bash
# Run split OTU
while read line; do

    VAL1="$( echo $line | cut -f1 -d " " )";
    VAL2="$( echo $line | cut -d " " -f2 )";
    Reformatre_trees.sh $VAL1 $VAL2;

done < burn
```

## Single Nuclotide Variant calling
Substitution calling is managed by a suite of Bash and Perl scripts, executed using the The_Run.sh package. This script must be executed in the same directory as the subs_selected.trees files, along with all required complementary scripts and files.

- combine_matrixes.sh
- Analyse_iter_ancestor_matrix6.sh
- Make_a_combined_matrix_universal.pl 
- extract_ancetre_seq_9.0.pl
- extract_ancetre_seq_12.0.pl
- Full_seq_length_line.txt
- Infos.txt
- Make_run.tt


```bash
# Run SNPs calling
The_Run.sh
```
## Figure representation

## Simulations 

## Non Bayesian Phylogenetic tree inference
