# HCMV Transcriptome Analysis Pipeline
**Author:** Peyton Fay  

## Project Overview
This pipeline provides an automated, end-to-end analysis of **Human Cytomegalovirus (HCMV)** transcriptomic data. The goal is to identify differentially expressed genes in human fibroblasts during infection and perform a *de novo* assembly of the viral genome.

### The Samples
The analysis uses RNA-Seq data from human fibroblasts infected with HCMV, obtained from the NCBI Sequence Read Archive (SRA).
* **Sample IDs:** `SRR5660030`, `SRR5660033`, `SRR5660044`, `SRR5660045`
* **Data Handling:** To ensure computational efficiency and meet performance requirements, raw FASTQ files were subset to a representative depth while maintaining biological signal for pipeline verification.

---

## 🔄 Analysis Workflow
The pipeline is managed by **Snakemake** and follows these logical steps:

1.  **Reference Construction:** Extracting 169 Coding Sequences (CDS) from the HCMV genome (NC_006273.2).
2.  **Quantification:** Building a **Kallisto** index and quantifying transcript abundance.
3.  **Differential Expression:** Running **Sleuth** in R to identify significant viral gene expression.
4.  **Enrichment/Filtering:** Using **Bowtie 2** to extract only viral reads, effectively removing host (human) DNA contamination.
5.  **Assembly:** Performing *de novo* assembly of filtered reads into contigs using **SPAdes**.
6.  **Verification:** Validating assembled contigs via **BLAST** against a Betaherpesvirinae database.



---

## 🛠 Technical Implementation (Terminal Commands)

### 1. Repository Management & Troubleshooting
During development, specific Git workflows were used to resolve merge conflicts and ensure repository integrity:
```bash
# Force-syncing local progress to GitHub to resolve conflicts
git add .
git commit -m "Finalized pipeline and documentation"
git push origin main --force
