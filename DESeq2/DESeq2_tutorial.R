###################################
#
#      Karl Grieshop
#      University of East Anglia
#      Data Science and Bioinformatics (BIO-7051B)
#      DESeq2 Tutorial - Differential Gene Expression Analysis
#      Data were prepared using MaleLimitedEvo_Pipeline.sh
#      2025-02-27
# 
###################################

# rm(list=ls()) # Clears the environment

# Packages
##########
# Function to check if a package is installed, and install it if not
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# Check and install required packages
install_if_missing("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
if (!requireNamespace("vsn", quietly = TRUE)) {
  BiocManager::install("vsn")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("clipr", quietly = TRUE)) {
  install.packages("clipr")
}

# Load necessary libraries
library(DESeq2)
library(vsn)
library(dplyr)
library(ggplot2)
library(clipr)

##########



# Now get chromosome locations from reference genome (will be useful later)
##########

all.genes <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/all.genes.tsv", sep="\t", header=FALSE)
colnames(all.genes) = c("FlyBaseID")

Xchr <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/X.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Xchr) = c("FlyBaseID")
Xchr$Chr <- rep("X", dim(Xchr)[1])

Ychr <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/Y.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Ychr) = c("FlyBaseID")
Ychr$Chr <- rep("Y", dim(Ychr)[1])

chr2L <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/2L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2L) = c("FlyBaseID")
chr2R <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/2R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2R) = c("FlyBaseID")
#
chr2 <- rbind(chr2L, chr2R)
chr2$Chr <- rep("2", dim(chr2)[1])

chr3L <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/3L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3L) = c("FlyBaseID")
chr3R <- read.delim(file="~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ExtractGenes/GeneLists/3R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3R) = c("FlyBaseID")
#
chr3 <- rbind(chr3L, chr3R)
chr3$Chr <- rep("3", dim(chr3)[1])

Chrs <- rbind(Xchr, Ychr, chr2, chr3) # Not all genes; just X, Y, 2, and 3.
Chrs$Chr <- as.factor(Chrs$Chr)
##########



# Set path to RNA-Seq read count data files
Data_path <- "~/Documents/BIO UEA/Teaching/Module - Data Science and Bioinformatics/Session10/RNA-Seq_Tutorial/ReadCounts/"
setwd(Data_path)

# List the files by category
files <- list.files(Data_path)

# Chop off the ".tsv" from those file names
samplename <- gsub('.{4}$', '', files)

# Set up factors for condition (must match order that files are listed; see "samplename" and "files")
sex <- factor(rep(c("Female", "Male"), each = 12))
geno <- factor(c(rep(c("NR", "Red"), each = 6, times = 2)))
rep <- factor(c(as.character(c(rep(1:6, each = 1, times = 4) ))))

# Create the sample table
sampleTable <- data.frame(sampleName = samplename, 
                          fileName = files, 
                          sex = sex,
                          geno = geno,
                          rep = rep)

# Display the structure of the sample table to check that factors are correctly set
str(sampleTable)

##########

# Subset data
A.f <- sampleTable[(sampleTable$sex == "Female"),] 
A.f$sex <- droplevels(A.f$sex)

A.m <- sampleTable[(sampleTable$sex == "Male"),]
A.m$sex <- droplevels(A.m$sex)

##########

# Set up contrast designs
dds.A.f.geno <- DESeqDataSetFromHTSeqCount(sampleTable = A.f, 
                                           design = ~ rep + geno) # model
dds.A.m.geno <- DESeqDataSetFromHTSeqCount(sampleTable = A.m, 
                                           design = ~ rep + geno) # model

# Set parameters for the focal contrast
minCountPerSample = 1 # you decide
minAvgPerCat = 10 # you decide
focal.contrast <- dds.A.f.geno # change accordingly (note confusing overwrite in DESeq2 documentation in this section)

# Specify the samples for each category of the focal contrasts
numerator <- samplename[7:12] # print this to make sure it's right! (e.g. Red females: samplename[7:12])
denominator <- samplename[1:6] # print this to make sure it's right! (e.g. NR females: samplename[1:6])

# Analysis details
factor.numerator.denominator = c("geno", "Red", "NR") # used later (change accordingly)
alpha.threshold = 0.05 # defaults to 0.05 but you could change it

# Set up count dataframe for focal contrast
countdf = DESeq2::counts(focal.contrast) # See Environment, "Large matrix", how large? Correct?)
# Filtering (set parameters above)
lowCountAnySample.numerator = sapply(1:(dim(countdf)[1]), function (x) prod(countdf[x, numerator]) < minCountPerSample)
lowCountAnySample.denominator = sapply(1:(dim(countdf)[1]), function (x) prod(countdf[x, denominator]) < minCountPerSample)
lowCountAnySample.Either =  lowCountAnySample.numerator | lowCountAnySample.denominator
avgCounts.numerator = rowMeans(countdf[, numerator])
avgCounts.denominator = rowMeans(countdf[, denominator])
goodAvgCount.numerator = avgCounts.numerator > minAvgPerCat
goodAvgCount.denominator = avgCounts.denominator > minAvgPerCat
goodAvgCount.Both = goodAvgCount.numerator & goodAvgCount.denominator
keep.these = goodAvgCount.Both & (!lowCountAnySample.Either)
# Filter the focal contrast dataframe 
focal.contrast.filtered = focal.contrast[keep.these]

# Do the analysis for that focal contrast using those filtered data
DESeq.Analysis = DESeq(focal.contrast.filtered)

# Get the results of that analysis
DESeq.Results = results(DESeq.Analysis, contrast = factor.numerator.denominator, alpha = alpha.threshold, independentFiltering=T)

# Remove Y genes (relevant to comparing males and females)
DESeq.Results$FlyBaseID = rownames(DESeq.Results)
DESeq.Results <- DESeq.Results[!(DESeq.Results$FlyBaseID %in% Ychr$FlyBaseID),]

# Remove genes on Chr 4 and pseudogenes (optional)
DESeq.Results$FlyBaseID = rownames(DESeq.Results)
DESeq.Results <- DESeq.Results[(DESeq.Results$FlyBaseID %in% Chrs$FlyBaseID),]

# Look at the metadata for DEseq's independent filtering 
plot(metadata(DESeq.Results)$filterNumRej, type="b", ylab="number of rejections", xlab="quantiles of filter")
lines(metadata(DESeq.Results)$lo.fit, col="red")
abline(v=metadata(DESeq.Results)$filterTheta)

# Have a basic look 
summary(DESeq.Results)
DESeq2::plotMA(DESeq.Results, ylim=c(-5,5), colSig = "red", colNonSig = "lightgray")

##########

# Perform PCA
vsd <- vst(focal.contrast, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("geno", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot PCA
ggplot(pcaData, aes(PC1, PC2, color=rep, shape=geno)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

##########

# Organize results into data frame
DESeq.Results$FlyBaseID = rownames(DESeq.Results)
Results.df <- data.frame(cbind(DESeq.Results$log2FoldChange, 
                               DESeq.Results$lfcSE, 
                               DESeq.Results$padj,
                               DESeq.Results$FlyBaseID))
colnames(Results.df) <- c("exp_geno", "se_geno", "padj", "FlyBaseID") # must be same order as previous line
Results.df$exp_geno <- as.numeric(Results.df$exp_geno)
Results.df$se_geno <- as.numeric(Results.df$se_geno)
Results.df$padj <- as.numeric(Results.df$padj)

##########

# Function to sort significant and non-significant genes by adding logical column
assign_sig <- function(contrast_df){
  contrast_df <- na.omit(contrast_df)
  contrast_df$Sig = FALSE
  for(i in 1:nrow(contrast_df)){
    if(contrast_df$padj[i] < alpha.threshold){
      contrast_df$Sig[i] = TRUE
    }
  }
  print(dim(contrast_df[contrast_df$Sig == TRUE, ])) # function will print...
  return(contrast_df)
}

# Assign significant genes
Results.df <- assign_sig(Results.df) # ... number of sig genes here ...
dim(Results.df[Results.df$Sig == TRUE, ]) # ... which should match this, correct?

##########

# Save the raw results (careful to name accordingly)
write.table(Results.df, file = "~/where/you/want/to/save/A.f.geno_raw.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)




##########
#
# Ok, so you've done the female Red/Non-red contrast, now do the males.
#
# Save that Results.df somewhere logical, named something logical. 
#
# To save a dataframe, or update/overwrite a saved dataframe, run:
# write.table(Results.df, file = "~/where/you/want/to/save/filename.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# Then...
# Go back up to line 91-95 (or thereabouts) and change:
# focal.contrast <- dds.A.f.geno 
# to 
# focal.contrast <- dds.A.m.geno 
# 
# AAAANNNNDDD what else needs changing in that area?
# 
#
# Once, coded for the male data, run back through the code to generate the male Red/Non-red contrast.
# 
# Again, save that male Red/Non-red contrast under a different name (don't overwrite the female one). 
#
##########



### QUESTION 1 ###
# 
# What do you notice that's different about the male Red/Non-red contrast?
#
##################


# To load a saved DESeq2 results dataframe, uncomment and run:
# Name.it <- read.delim(~/where/you/saved/it/filename.tsv", sep = '\t', header = TRUE))
# So, for example, if I called those "A.f.geno_raw.tsv" and "A.f.geno_raw.tsv" maybe:
A.f.geno <- read.delim("where/you/saved/A.f.geno_raw.tsv") 
A.m.geno <- read.delim("where/you/saved/A.m.geno_raw.tsv")



##########
#
# Ok, now lets focus on the candidate "differentially expressed" genes.
#
# How could you isolate those in females?
#
# What about males?
#
# What if you want to compare the female and male candidates?
#
##########



# For females, just:
A.f.geno.can <- A.f.geno[A.f.geno$Sig == "TRUE",]
# Save to file
write.table(A.f.geno.can, file = "~/where/you/want/to/save/A.f.geno_candidates.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

# But if we do that for males... 
A.m.geno.can <- A.m.geno[A.m.geno$Sig == "TRUE",] # Unsatisfying :-/ 
# Don't save that.



### QUESTION 2 ###
# 
# What's your intuition: Can we compare the characteristics of a set of 353 genes to that of 16 genes?
#
##################




### QUESTION 3 ###
# 
#  What should you do in this situation? Discuss in groups.
#
##################




# One solution:
# This is only for the Exp. males comparison. Because the 5% FDR threshold is too stringent (due to variation between replicate measures), 
# we take 350 genes with the highest log2FC value instead (since that's how many we have for females).

# Be sure that "Results.df" is what you think it is (should be the male constrast results for this)
Results.df <- Results.df[order(abs(Results.df$exp_geno), decreasing = T),] # order from greatest log2FC value
colnames(Results.df)[colnames(Results.df) == "Sig"] = "Top.Sig" # keep the 5% FDR assigned significant genes, but assign it another name
# Re-assign the significance column, and initially set to FALSE
Results.df$Sig <- FALSE # change the 350 with greates log2FC to $Sig == "TRUE"
i = 1 # start count for number of significant gene
# go through the list of genes until 350 genes are assigned as significant in the male data
while(dim(Results.df[Results.df$Sig,])[1] < 350){
  if(!Results.df$FlyBaseID[i] %in% Ychr$geneID & # exclude Y-linked genes
     Results.df$FlyBaseID[i] %in% Chrs$FlyBaseID){ # exclude genes on Chr 4 and pseudogenes
  Results.df$Sig[i] <- TRUE # changes the 350 with greatest log2FC to $Sig == "TRUE"
  }
  i = i + 1 # update count
}
droplevels(Results.df)
dim(Results.df[Results.df$Sig,]) # ensure there is 350 here.

# Take just those top 350
A.m.geno.can <- Results.df[Results.df$Sig,]
# Save it
write.table(A.m.geno.can, file = "~/where/you/want/to/save/A.m.geno_candidates.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

##########


# If you haven't already, load raw results 
A.f.geno <- read.delim("where/you/saved/A.f.geno_raw.tsv")
A.m.geno <- read.delim("where/you/saved/A.m.geno_raw.tsv")
# Combine results for males and females populations (all genes, background set)
Exp.geno.raw <- merge(A.m.geno, A.f.geno, by = "FlyBaseID", all = TRUE)


# If you haven't already, load candidate results
A.f.geno.can <- read.delim("where/you/saved/A.f.geno_candidates.tsv")
A.m.geno.can <- read.delim("where/you/saved/A.m.geno_candidates.tsv")

# Combine the candidates for males and females (removing NAs)
Exp.geno.can <- merge(A.m.geno.can, A.f.geno.can, by = "FlyBaseID", all = TRUE)
colnames(Exp.geno.can) <- c("FlyBaseID", "A.m.exp_geno", "A.m.se_geno", "A.m.padj", "A.m.TopSig", "A.m.Sig", # Check order of these is correct
                            "A.f.exp_geno", "A.f.se_geno", "A.f.padj", "A.f.Sig")
Exp.geno.can <- Exp.geno.can %>% mutate(Sig = ifelse(!is.na(A.m.Sig) & A.m.Sig, TRUE,
                                                     ifelse(!is.na(A.f.Sig) & A.f.Sig, TRUE, 
                                                            ifelse(is.na(A.m.Sig) & is.na(A.f.Sig), NA, FALSE)))) 
Exp.geno.can <- Exp.geno.can[!is.na(Exp.geno.can$Sig),]

# Save to file
write.table(Exp.geno.can, file = "~/where/you/want/to/save/All.geno_candidates.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)



# Combine the background for males and females (removing NAs)
Exp.geno.background <- merge(A.m.geno, A.f.geno, by = "FlyBaseID", all = FALSE)
Exp.geno.background <- Exp.geno.background[!is.na(Exp.geno.background$Sig),]



### QUESTION 4 ###
# 
# Of the differentially expressed genes in males and females, only 40 overlap. 
# Did male-limited selection cause a correlated change in female gene expression?
# How did authors test this?
# Hint: got to paper, ctl+f, "overlap"
#
##################




##########
#
# GO enrichment for candidate differentially expressed genes
# 
##########

# Target set
write_clip(Exp.geno.can$FlyBaseID) 
# Background set
write_clip(Exp.geno.background$FlyBaseID) 
# Background set
# These were used here: https://biit.cs.ut.ee/gprofiler/gost
# Alternatively, use: https://bioinformatics.sdstate.edu/go/
##########




### QUESTION 5 ###
# 
# Could you think of an alternative GO enrichment analysis on a different
# set of candidate genes? Try it.
#
##################


