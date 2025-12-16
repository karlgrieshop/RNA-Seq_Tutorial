# RNA-Seq_Tutorial for Data Science and Bioinformatics (BIO-7051B) Session 5
STAR | HTSeq | DESeq2 workflow 

---

## Learning objectives
- Understand the pipeline steps to go from FastQ → read counts (STAR | HTSeq).
- Run basic DESeq2 workflows in R and interpret simple differential expression results.

---

## Origin
This tutorial is modified from select subdirectories of the public repository [MaleLimitedEvo](https://github.com/mchlleliu/MaleLimitedEvo), which accompanies the publication Grieshop et al. (2025) in *Molecular Biology and Evolution (MBE)*. It contains scripts and pipelines for processing RNA-Seq data, including extracting gene lists, downloading reference genomes, renaming files, and generating read-count TSV files for downstream analysis using DESeq2 in R.

---

## Background
In six replicate populations of 1000 flies, a dominant marker (DsRed) on Chromosome 2 was used to force a “Red” pool of genetically variable chromosomes through exclusive father-to-son inheritance, while a complimentary pool of “NonRed” chromosomes was inherited primarily from mothers to daughters. After 100 generations, males carrying a Red chromosome copy exhibited greater mating success than males with only NonRed chromosomes, consistent with the accumulation of male-benefit/female-detriment sexually antagonistic alleles in the Red pool relative to NonRed. We analysed differentially expressed genes between flies with and without Red chromosomes.

---

## Two paths: Step 1&2, or just Step 2

**Step 1 (optional, advanced, begin *prior* to Session5)**
- Try to run the full /Pipeline 
- Requires all elements of this repo (*listed in order of execution under `Contents`, below*).
- A recommended workflow:
  - (As in BED_overlap) `fork` the repo → `clone` it locally → Edit files → `push` changes → `clone` or `pull` to HPC workspace → run scripts.
  - Alternatives: 
   - `scp` essential local changes to HPC
   - Work entirely on the HPC: `clone` to HPC → Edit directly in HCP (using `vim` or `nano`).
- Once you have TSV read-count files, move onto Option 2 (next).

**Step 2 (required, beginner, begin *during* Session5)**
- If you skipped Step 1, you can `scp` the TSV read-count files from the HPC to your local workspace (see `ReadCounts/`), and 
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

### 7. `Plotting/`
- **Script**: `DE_Plotting_data.R`, `DE_SBGEdist_Fig2.S8.S9.R`, `Mishra_et.al_SBGE.R`
- **Data**: `DifferentialGeneExpression.whole.bodies.tsv`
- **Purpose**: *NOT* a core component of the tutorial. Starting material for Group Project idea. Incomplete code and data for visualising results.
- **Execution**:  data wrangling and ggplot in R, using outputs from DESeq2 results table.
- **Output**: Results figures.

---

## Group Project ideas

Take this further in some way or ways (remember: chief marking criteria for group project is the reproducibility):
- **Pipeline work:** Convert the FastQProcessing/Pipeline to a [Nextflow pipeline](https://www.nextflow.io/rna-seq-pipeline.html) or a [Snakemake pipeline](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/rna-seq-star-deseq2.html).
  - Take it further, compare outputs (e.g. are read-count TSV files correlated to those of FastQProcessing/Pipeline?).
  - Take it further, compare analysis (e.g. do different pipelines result in different expression results?).
- **Outliers:** Remove outlier sample. Look at PCA among samples in R. Do any samples stand out? If so, try to exclude that/those sample(s) from the analysis and re-run. Compare results. Did they change?
- **Data Visualisation:** `Plotting/` contains incomplete starting material to get you started with some data visualisation ideas for how to represent the results. See plots from Grieshop *et al*. 2025 for inspiration, or to see what they haven't done. 
- **Statistical:** Analyse the chromosomal distribution of differentially expressed genes. Are there more on X, 2 or 3 than expected by chance? Does this make sense? Recall the focal manipulated chromsome of the experiment was Chromosome 2.  

---

## Contact & Questions
For questions about this tutorial, contact:
Karl Grieshop  
School of Biological Sciences  
University of East Anglia  
k.grieshop@uea.ac.uk