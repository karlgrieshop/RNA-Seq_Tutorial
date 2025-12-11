# RNA-Seq_Tutorial/DESeq2/README.md

## This is local work - Open the R script locally and follow it
It's best and easiest to run R Studio locally, and the .tsv files are small enough to download to your local machine. 

### If you don't have the TSV read-count files locally already, look at RNA-Seq_Tutorial/ReadCounts/README.md to see where to find the TSV files. Then, use scp (secure shell copy) to get the files on your local machine, like so:
- Open a local terminal session
```bash
cd /<where_you_want_them/>
# `scp` is secure copy protocol, and allows you to `cp` between HPC and local workspaces:
scp abc12xyz@hali.uea.ac.uk:~/<where_you_have_them/*.tsv <local_directory_name>
# And enter your UEA HPC password when prompted,
# Wait for download to complete before closing shell.
```

### Options for getting the DESeq2_tutorial.R script to your local machine
- You may have already have (it's in the repo, if you cloned it locally). 
- Or you could download the specific R sript from the GitHub RNA-Seq_Tutorial/.
- Or you could use `scp` as shown above.

