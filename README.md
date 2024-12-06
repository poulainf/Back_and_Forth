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
These commands will produce a file named Human_Viruses_RefSeq_and_neighbors_genome_data-cleaned.tsv. This file associates each full genome sequence of a human virus from the NCBI nucleotide database with a specific reference sequence.

The cleaned table file can then be used to download FASTA and GenBank files using esearch from the Entrez Direct Utilities.
