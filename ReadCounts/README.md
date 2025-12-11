# RNA-Seq_Tutorial/ReadCounts/README.md

**There are not .tsv read-count files stored in this GitHub repo.** 

**If you are a student (or otherwise) in Data Science and Bioinformatics (BIO-7051B),** you can copy the data from the shared HPC directory /gpfs/data/BIO-DSB/Session5/ to your ~/scratch/ directory or some other logically named location off of your ~/ directory, e.g.:

```bash
ssh <abc12xyz>@hali.uea.ac.uk
interactive-bio-ds
cp /gpfs/data/BIO-DSB/Session5/RNA-Seq_Tutorial/ReadCounts/*.tsv ~/scratch/<where_you_wish>
```

**If you do not have access to /gpfs/data/BIO-DSB/Session5/** you can apply the code in Pipeline/ to the relevant original fastq files publicly available on the [sequence read archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1184789). **Note** that this tutorial excludes the Control samples in the SRA from the Grieshop et al. (2025) (PRJNA1184789) (job array is for 24 files, not 36).