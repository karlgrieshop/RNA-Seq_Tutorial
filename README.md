# RNA-Seq_Tutorial
STAR | HTSeq | DESeq2 workflow for transcriptomics tutorial for Data Science and Bioinformatics (UEA, BIO7051B)

## Origin
This tutorial is modified from select subdirectories of the  public repository [MaleLimitedEvo](https://github.com/mchlleliu/MaleLimitedEvo), which accompanies the publication Grieshop et al. (2025) in *Molecular Biology and Evolution (MBE)*. It contains scripts and pipelines for processing RNA-Seq data, including extracting gene lists, downloading reference genomes, renaming files, and generating read-count TSV files for downstream analysis using DESeq2 in R. 

## Background
In six replicate populations of 1000 flies, a dominant marker (DsRed) on Chromosome 2 was used to force a “Red” pool of genetically variable chromosomes through exclusive father-to-son inheritance, while a complimentary pool of “NonRed” chromosomes was inherited primarily from mothers to daughters. After 100 generations, males carrying a Red chromosome copy exhibited greater mating success than males with only NonRed chromosomes, consistent with the accumulation of male-benefit/female-detriment sexually antagonistic alleles in the Red pool relative to NonRed. We analysed differentially expressed genes between flies with and without Red chromosomes. 

## Tutorial overview
This directory contains materials for Session 10 of the Data Science and Bioinformatics module at the University of East Anglia. The focus of this session is to:
1. Document how RNA-Seq FastQ files were processed to read count .tsv files.
2. Demonstrate how to analyse read counts using DESeq2 in R.

## Contents

### 1. `Pipeline/`
- **Script**: `FastQ-to-ReadCounts_Pipeline.txt`
- **Purpose**: Processes paired-end RNA-Seq data from FastQ files into read-count TSV files. This includes genome indexing, read alignment, BAM file processing, and read counting.
- **Execution**: Designed for interactive execution locally or on a server or an HPC.
- **Dependency**: Requires the reference genome and annotation files generated by the `GetReference/Get_Reference_Genome.sh` script.
- **Output**: Read-count TSV files for each sample.
- **Notes**: See `fastqFileRenaming/`

### 2. `ExtractGenes/`
- **Script**: `Extract_Genes.sh`
- **Purpose**: Extracts gene lists for specific chromosomes (e.g., X, Y, autosomes, mitochondrial genome) from a GTF file, used for filtering and characterization in MaleLimitedEvo/DESeq2Run
/DESeq2Run.R.
- **Execution**: Runs on an HPC using SLURM.
- **Output**: Chromosome-specific gene lists in TSV format.

### 3. `GetReference/`
- **Script**: `Get_Reference_Genome.sh`
- **Purpose**: Downloads and prepares the Drosophila melanogaster reference genome and annotation files. It also generates indices for BWA, samtools, and STAR.
- **Execution**: Runs on an HPC using SLURM.
- **Output**: Reference genome files, GTF annotation file, and genome indices.

### 4. `tsvFileRenaming/`
- **Script**: `renameTSVfiles.sh`
- **Purpose**: Renames TSV files based on their sample name components for easier sorting and analysis in R.
- **Execution**: Runs on an HPC using SLURM.
- **Output**: Renamed TSV files.

### 5. `fastqFileRenaming/`
- **Script**: `renameFASTQfiles.sh`
- **Purpose**: Renames FastQ files based on a predefined mapping for consistency and clarity.
- **Execution**: Does not require SLURM; can be run locally or on a server or an HPC.
- **Output**: Renamed FastQ files.

### 6. `DESeq2/`
- **Script**: `DESeq2_tutorial.R`
- **Purpose**: Contains the DESeq2 tutorial script for differential gene expression analysis.
- **Execution**: DESeq2 in R, using the data in `ReadCounts/` prepared using the `Pipeline/FastQ-to-ReadCounts_Pipeline.txt` commands.
- **Output**: Results of differential gene expression analysis.

### 7. `ReadCounts/`
- **Purpose**: Contains the read count files generated from `Pipeline/FastQ-to-ReadCounts_Pipeline.txt` for the analysis documented in `DESeq2/DESeq2_tutorial.R`.
- **Notes**: See `tsvFileRenaming/`.

## Software Requirements
- **STAR**: v2.7.9a or later
- **GATK**: v4.2.3.0-1 or later
- **HTSeq**: v0.13.5 or later
- **GNU parallel**: Required for parallelizing sample processing
- **SLURM**: Required for job scheduling on an HPC
- **rename**: Required for renaming files in `renameFASTQfiles.sh`

## Reference Genome and Annotation
- **Reference Genome**: Drosophila melanogaster BDGP6.28, release 102
- **Annotation File**: GTF file from the same release

## Notes
- Ensure all required software is installed and accessible in your environment before running the scripts.
- Some scripts are designed to run on an HPC using SLURM. Modify the SLURM directives as needed for your system.
- The `Pipeline/` script was executed interactively on a private server but can be adapted for HPC use.

## Contact
For questions or issues, contact [Karl Grieshop](https://github.com/karlgrieshop).
