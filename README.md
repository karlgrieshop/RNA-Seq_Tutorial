# RNA-Seq_Tutorial
STAR | HTSeq | DESeq2 workflow for transcriptomics tutorial for Data Science and Bioinformatics (UEA, BIO7051B)

## PLEASE 

** Do not share this data outside this module.**

## Origin
The data stem from the unpublished article:
Grieshop, K., Liu, M. J., Frost, R. S., Lindsay, M. P., Bayoumi, M., Brengdahl, M. I., Molnar, R. I. & Agrawal, A. F. (2024). Expression divergence in response to sex-biased selection. bioRxiv, 2024-11.

## Background
In six replicate populations of 1000 flies, a dominant marker (DsRed) on Chromosome 2 was used to force a “Red” pool of genetically variable chromosomes through exclusive father-to-son inheritance, while a complimentary pool of “NonRed” chromosomes was inherited primarily from mothers to daughters. After 100 generations, males carrying a Red chromosome copy exhibited greater mating success than males with only NonRed chromosomes, consistent with the accumulation of male-benefit/female-detriment sexually antagonistic alleles in the Red pool relative to NonRed. We analysed differentially expressed genes between flies with and without Red chromosomes. 

## Tutorial overview
This directory contains materials for Session 10 of the Data Science and Bioinformatics module at the University of East Anglia. The focus of this session is to:
1. Document how RNA-Seq FastQ files were processed to read count .tsv files.
2. Demonstrate how to analyse read counts using DESeq2 in R.

## Contents
- `DESeq2/`: Contains the DESeq2 tutorial script for differential gene expression analysis.
- `ExtractGenes/`: Contains the script for extracting gene lists from the GTF file.
- `FastQProcessing/`: Contains subdirectories and scripts for processing FastQ files.
  - `GetReference/`: Scripts for downloading and indexing the reference genome.
  - `Pipeline/`: Scripts for processing RNA-Seq data from FastQ files to .tsv files.
  - `RenamingFiles/`: Scripts for renaming .tsv files for easier sorting and analysis in R.
- `ReadCounts/`: Contains the read count files generated from the RNA-Seq data.

## Detailed Descriptions

### DESeq2
This directory contains the `DESeq2_tutorial.R` script, which demonstrates how to perform differential gene expression analysis using DESeq2 in R. The data were prepared using the `MaleLimitedEvo_Pipeline.sh` script.

### ExtractGenes
This directory contains the `Extract_Genes.sh` script, which extracts gene lists from the GTF file for different chromosomes.

### FastQProcessing
This directory contains subdirectories and scripts for processing FastQ files.

#### GetReference
Contains the `Get_Reference_Genome.sh` script, which downloads and indexes the reference genome files required for RNA-Seq analysis.

#### Pipeline
Contains the `MaleLimitedEvo_Pipeline.sh` script, which processes RNA-Seq data from FastQ files to .tsv files using the indexed reference genome.

#### RenamingFiles
Contains the `renameFiles.sh` script, which renames `.tsv` files based on specific components extracted from their original filenames for easier sorting and analysis in R.

### ReadCounts
This directory contains the resultant read count files generated from the RNA-Seq data.

## Usage

### Running the Scripts

1. **Download and Index the Reference Genome:**
   ```bash
   sbatch FastQProcessing/GetReference/Get_Reference_Genome.sh
   ```

2. **Process RNA-Seq Data:**
   ```bash
   sbatch FastQProcessing/Pipeline/MaleLimitedEvo_Pipeline.sh
   ```

3. **Rename .tsv Files:**
   ```bash
   sbatch FastQProcessing/RenamingFiles/renameFiles.sh
   ```

4. **Extract Gene Lists:**
   ```bash
   sbatch ExtractGenes/Extract_Genes.sh
   ```

5. **Perform Differential Gene Expression Analysis:**
   Open and run the `DESeq2/DESeq2_tutorial.R` script in R.

## Notes
- Ensure that the input directories and file paths are correctly specified in the scripts before running them.
- The `MaleLimitedEvo_Pipeline.sh` script has not been tested yet and may require adjustments.

## Contact
For any issues or questions, please contact K.Grieshop@uea.ac.uk