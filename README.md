# RNA-Seq_Tutorial
STAR | HTSeq | DESeq2 workflow for transcriptomics tutorial for Data Science and Bioinformatics (UEA, BIO-7051B)

---

## Origin
This tutorial is modified from select subdirectories of the public repository [MaleLimitedEvo](https://github.com/mchlleliu/MaleLimitedEvo), which accompanies the publication Grieshop et al. (2025) in *Molecular Biology and Evolution (MBE)*. It contains scripts and pipelines for processing RNA-Seq data, including extracting gene lists, downloading reference genomes, renaming files, and generating read-count TSV files for downstream analysis using DESeq2 in R.

---

## Background
In six replicate populations of 1000 flies, a dominant marker (DsRed) on Chromosome 2 was used to force a “Red” pool of genetically variable chromosomes through exclusive father-to-son inheritance, while a complimentary pool of “NonRed” chromosomes was inherited primarily from mothers to daughters. After 100 generations, males carrying a Red chromosome copy exhibited greater mating success than males with only NonRed chromosomes, consistent with the accumulation of male-benefit/female-detriment sexually antagonistic alleles in the Red pool relative to NonRed. We analysed differentially expressed genes between flies with and without Red chromosomes.

---

## Learning objectives
- Understand the pipeline steps to go from FastQ → read counts (STAR | HTSeq).
- Run basic DESeq2 workflows in R and interpret simple differential expression results.

---

## Steps

**Step 1 (optional, advanced, HPC, begin *prior* to Session5)**
- Try to run the full /Pipeline 
- Requires all elements of this repo (*listed in order of execution under `Contents`, below*).
- A recommended workflow:
  - (As in BED_overlap) `fork` the repo → `clone` it locally → Edit files → `push` changes → `clone` or `pull` to HPC workspace → run scripts.
  - Alternatives: 
   - `scp` essential local changes to HPC
   - Work entirely on the HPC: `clone` to HPC → Edit directly in HCP (using `vim` or `nano`).
- Once you have TSV read-count files, move onto Option 2 (next).

**Step 2 (required, local, R, complete *during* Session5)**
- If you skipped Step 1, you can `scp` the TSV read-count files from the HPC to your local workspace (see `ReadCounts/`).
- Analyse read-counts in R locally `DESeq2/`.

---

## Contents

### 1. `FastQProcessing/GetReference/`
- **Script**: `Get_Reference_Genome.sh`
- **Purpose**: Downloads and prepares the *Drosophila melanogaster* reference genome and annotation files. It also generates indices for BWA, samtools, and STAR.
- **Execution**: Runs on an HPC using SLURM.
- **Output**: Reference genome files, GTF annotation file, and genome indices.

### 2. `FastQProcessing/ExtractGenes/`
- **Script**: `Extract_Genes.sh`
- **Purpose**: Extracts gene lists for specific chromosomes (e.g., X, Y, autosomes, mitochondrial genome) from a GTF file, used for filtering and characterization in DESeq2 analysis.
- **Execution**: Runs on an HPC using SLURM.
- **Output**: Chromosome-specific gene lists in TSV format.

### 3. `FastQProcessing/Pipeline/`
- **Scripts**:  
  - `FastQ-to-ReadCounts_Pipeline.sh` 
- **Purpose**: Processes paired-end RNA-Seq data from FastQ files into read-count TSV files (STAR | BWA | samtools | HTSeq). 
- **Execution**:  
  - `sbatch` FastQ-to-ReadCounts_Pipeline.sh
- **Dependency**: 
  - Reference genome and annotation files (from `FastQProcessing/GetReference/).
  - Conda environment (see notes in `FastQ-to-ReadCounts_Pipeline.sh`).
- **Output**: Read-count TSV files for each sample.

### 4. `FastQProcessing/tsvFileRenaming/`
- **Script**: `renameTSVfiles.sh`
- **Purpose**: Renames TSV files (dangerous) in a careful, reproducible way, for easier sorting and analysis in R.
- **Execution**: Runs on an HPC using SLURM.
- **Output**: Renamed TSV files.

### 5. `ReadCounts/`
- **Purpose**: Contains the read count files generated from `Pipeline/` for the analysis documented in `DESeq2/DESeq2_tutorial.R`.
- **Notes**: See `tsvFileRenaming/`.

### 6. `DESeq2/`
- **Script**: `DESeq2_tutorial.R`
- **Purpose**: Contains the DESeq2 tutorial script for differential gene expression analysis.
- **Execution**: DESeq2 in R, using the data in `ReadCounts/` prepared using the `Pipeline/` scripts.
- **Output**: Results of differential gene expression analysis.
- **Notes**: You will run DESeq2 locally in R Studio, so you will need to copy the .tsv files somewhere locally and import them into R.

## Contact & Questions
For questions about this tutorial, contact:
Karl Grieshop  
School of Biological Sciences  
University of East Anglia  
k.grieshop@uea.ac.uk