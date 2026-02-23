#!/bin/bash

# Fay's HCMV Pipeline Execution Script
#this was done to avoid the issue of snakemake not being able to find the blast database when running the pipeline
echo "Starting HCMV Pipeline..."

# 1. Run the Snakemake pipeline
# Note: Ensure you have Snakemake installed and the Snakefile is in the current directory
#this terminal command will run the snakemake pipeline, which will execute all the steps defined in the Snakefile, including quality control, filtering, assembly, and BLAST identification
echo "Running Snakemake workflow..."
snakemake --cores 1 --latency-wait 60

# 2. Generate Read Counts (Step 4)
#this will count the number of reads remaining after filtering for each sample and save the counts to a text file
#i did this by counting the number of lines in the filtered FASTQ files (which contain the reads that passed the quality control and filtering steps) and dividing by 4 (since each read is represented by 4 lines in a FASTQ file)
echo "Generating post-filtering read counts..."
for s in SRR5660030 SRR5660033 SRR5660044 SRR5660045; do
    echo -n "$s After: " && expr $(cat results/${s}_mapped.1.fastq | wc -l) / 4
done > results/read_counts.txt

# 3. Build BLAST Database
#this will create a local BLAST database from the HCMV CDS sequences, which will be used for the BLAST identification step to identify the assembled contigs
echo "Building local BLAST database..."
makeblastdb -in hcmv_cds.fasta -dbtype nucl -out betaherpesvirinae

# 4. Run BLAST on the assembly (Step 6)
#this will take the assembled contigs from the assembly step and run a BLAST search against the local BLAST database created in the previous step to
#  identify the contigs based on sequence similarity. The results will be saved in a tabular format with specific fields such as subject accession, 
# percent identity, alignment length, etc.
echo "Running BLAST identification..."
blastn -query results/SRR5660033_assembly/contigs.fasta -db betaherpesvirinae -max_target_seqs 1 -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" > results/blast_results.txt

echo "Pipeline complete. Results are in the results/ folder."