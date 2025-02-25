#!/bin/bash
#SBATCH --time=1-0
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Extract_Genes
#SBATCH -o /gpfs/home/rpr23sxu/Teaching/MaleLimitedEvo/ExtractGenes/Output_Messages/Extract_Genes.out
#SBATCH -e /gpfs/home/rpr23sxu/Teaching/MaleLimitedEvo/ExtractGenes/Error_Messages/Extract_Genes.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rpr23sxu@uea.ac.uk

# Define the GTF file path
GTF_FILE="/gpfs/home/rpr23sxu/Teaching/MaleLimitedEvo/GetReference/References/Drosophila_melanogaster.BDGP6.28.102.gtf"

# X genes
awk '{if ($1 == "X") print $0;}' $GTF_FILE | awk '{print $10}' | sed 's/"//g' | sed 's/;//g' | uniq > X.chromosome.genes.tsv

# Y genes
awk '{if ($1 == "Y") print $0;}' $GTF_FILE | awk '{print $10}' | sed 's/"//g' | sed 's/;//g' | uniq > Y.chromosome.genes.tsv

# 2L genes
awk '{if ($1 == "2L") print $0;}' $GTF_FILE | awk '{print $10}' | sed 's/"//g' | sed 's/;//g' | uniq > 2L.chromosome.genes.tsv

# 2R genes
awk '{if ($1 == "2R") print $0;}' $GTF_FILE | awk '{print $10}' | sed 's/"//g' | sed 's/;//g' | uniq > 2R.chromosome.genes.tsv

# 3L genes
awk '{if ($1 == "3L") print $0;}' $GTF_FILE | awk '{print $10}' | sed 's/"//g' | sed 's/;//g' | uniq > 3L.chromosome.genes.tsv

# 3R genes
awk '{if ($1 == "3R") print $0;}' $GTF_FILE | awk '{print $10}' | sed 's/"//g' | sed 's/;//g' | uniq > 3R.chromosome.genes.tsv

# Mitochondrial genes
awk '{if ($1 == "mitochondrion genome") print $0;}' $GTF_FILE | awk '{print $10}' | sed 's/"//g' | sed 's/;//g' | uniq > mitochondrial.genes.tsv

# All genes (autosomes plus sex chromosomes)
awk '{if ($1 == "X" || $1 == "Y" || $1 == "3R" || $1 == "3L" || $1 == "2R" || $1 == "2L" || $1 == "4") print $0;}' $GTF_FILE | awk '{print $10}' | sed 's/"//g' | sed 's/;//g' | uniq > all.genes.tsv
