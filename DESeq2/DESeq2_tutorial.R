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

rm(list=ls()) # Clears the environment

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

# Load necessary libraries
library(DESeq2)
library(vsn)
library(dplyr)

##########

# Set path to data files
Data_path <- "~/Documents/Omics/DsRed/Transcriptomics/Subset/"
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
numerator <- samplename[7:12] # print this to make sure it's right! (e.g. simple: samplename[7:12]; e.g. specific samplename[c(7:12, 19:24)]  )
denominator <- samplename[1:6] # print this to make sure it's right! (e.g. simple: samplename[1:6]; e.g. specific samplename[c(1:6, 13:18)]  )

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

# Remove Y genes and DsRed marker genes
Ychr <- read.delim("~/Desktop/UofT/PRJ1/data/Y.chromosome.genes.tsv", sep = '\t', header = TRUE)
DsRed_genes <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/dmel_2R_DsRed_ids.tsv", header=FALSE)
DESeq.Results$FlyBaseID = rownames(DESeq.Results)
DESeq.Results <- DESeq.Results[!(DESeq.Results$FlyBaseID %in% Ychr$geneID),]
DESeq.Results <- DESeq.Results[!DESeq.Results$FlyBaseID %in% DsRed_genes$V1,]

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
colnames(Results.df) <- c("exp_geno", "se_geno", "padj", "FlyBaseID")
Results.df$exp_geno <- as.numeric(Results.df$exp_geno)
Results.df$se_geno <- as.numeric(Results.df$se_geno)
Results.df$padj <- as.numeric(Results.df$padj)

# Save to file
write.table(Results.df, file = "~/Desktop/UofT/SSAV_RNA/Results/A.m.geno_raw.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

##########

# Load DESeq2 result for comparison between Red and NonRed A-females
A.f.geno <- read.delim("Results/A.f.geno_raw.tsv")

# Function to sort significant and non-significant genes by adding logical column
assign_sig <- function(contrast_df){
  contrast_df <- na.omit(contrast_df)
  contrast_df$Sig = FALSE
  for(i in 1:nrow(contrast_df)){
    if(contrast_df$padj[i] < alpha.threshold){
      contrast_df$Sig[i] = TRUE
    }
  }
  print(dim(contrast_df[contrast_df$Sig == TRUE, ]))
  return(contrast_df)
}

# Assign significant genes
Results.df <- assign_sig(Results.df)
dim(Results.df[Results.df$Sig == TRUE, ])

# Save to file
write.table(Results.df, file = "~/Desktop/UofT/SSAV_RNA/Results/A.f.geno_candidates.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

##########

# Load DESeq2 result for comparison between Red and NonRed A-males
A.m.geno <- read.delim("Results/A.m.geno_raw.tsv")

# Function to sort significant and non-significant genes by adding logical column
assign_sig_male <- function(contrast_df){
  contrast_df <- na.omit(contrast_df)
  contrast_df$Sig = FALSE
  for(i in 1:nrow(contrast_df)){
    if(contrast_df$padj[i] < alpha.threshold){
      contrast_df$Sig[i] = TRUE
    }
  }
  print(dim(contrast_df[contrast_df$Sig == TRUE, ]))
  return(contrast_df)
}

# Assign significant genes
Results.df <- assign_sig_male(Results.df)
dim(Results.df[Results.df$Sig == TRUE, ])

# this is only for the Exp. males comparison. Because the 5% FDR threshold is too stringent (due to variation between replicate measures), 
# we take 350 genes with the highest log2FC value instead.
Results.df <- Results.df[order(abs(Results.df$exp_geno), decreasing = T),] # order from greatest log2FC value
colnames(Results.df)[colnames(Results.df) == "Sig"] = "Top.Sig" # keep the 5% FDR assigned significant genes, but assign it another name
# Re-assign the significance column, and initially set to FALSE
Results.df$Sig <- FALSE
i = 1 # start count for number of significant gene
# go through the list of genes until 350 genes are assigned as significant in the male data
while(dim(Results.df[Results.df$Sig,])[1] < 350){
  if(!Results.df$FlyBaseID[i] %in% Ychr$geneID & # exclude Y-linked genes
     Results.df$FlyBaseID[i] %in% Chrs$FlyBaseID){ # exclude genes on Chr 4 and pseudogenes
  Results.df$Sig[i] <- TRUE
  }
  i = i + 1 # update count
}
droplevels(Results.df)
dim(Results.df[Results.df$Sig,]) # ensure there is 350 here.

# Save to file
write.table(Results.df, file = "~/Desktop/UofT/SSAV_RNA/Results/A.m.geno_candidates.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

##########

# Load results
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")

# Combine results for experimental populations
Exp.geno <- merge(A.m.geno, A.f.geno, by = "FlyBaseID", all = TRUE)
colnames(Exp.geno) <- c("FlyBaseID", "A.m.exp_geno", "A.m.se_geno", "A.m.padj", "A.m.TopSig", "A.m.Sig",
                          "A.f.exp_geno", "A.f.se_geno", "A.f.padj", "A.f.Sig")
Exp.geno <- Exp.geno %>% mutate(Sig = ifelse(!is.na(A.m.Sig) & A.m.Sig, TRUE,
                                               ifelse(!is.na(A.f.Sig) & A.f.Sig, TRUE, 
                                                      ifelse(is.na(A.m.Sig) & is.na(A.f.Sig), NA, FALSE)))) 
Exp.geno <- Exp.geno[!is.na(Exp.geno$Sig),]

# Only keep concordant changes
Exp.geno.con <- na.omit(Exp.geno[(Exp.geno$A.f.exp_geno > 0 & Exp.geno$A.m.exp_geno > 0) |
                             (Exp.geno$A.f.exp_geno < 0 & Exp.geno$A.m.exp_geno < 0),])

# Save to file
write.table(Exp.geno, file = "~/Desktop/UofT/SSAV_RNA/Results/All.geno_candidates.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
