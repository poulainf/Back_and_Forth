# Back_and_Forth

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
The cleaned table file can then be used to download FASTA and GenBank files using `esearch` from the Entrez [Direct utilitiez](https://www.ncbi.nlm.nih.gov/books/NBK179288/).The following command allow an automated files download nad concatenated them. Sequence genomes clustering into Operative Taxonomic Units are next performed based on taxonomic information provided by GeneBank datas and Neighbor-to-RefSeq Sequence File. 

```bash
# Download the selected genome data files and GeneBank files from NCBI
Download_Fasta_GeneBank.pl Human_Viruses_RefSeq_and_neighbors_genome_data-cleaned.tsv

# Concatenate files
grep human Viruses_RefSeq_and_neighbors_genome_data.tab > Human_Viruses_RefSeq_and_neighbors_genome_data.tab

# Cluster sequence into OTU 
OTU_fasta_concatenation.pl Human_Viruses_RefSeq_and_neighbors_genome_data.tab
```
