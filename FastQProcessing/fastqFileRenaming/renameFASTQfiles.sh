#!/bin/bash

# This script renames FastQ files based on a predefined mapping.
# The order of the 'rename' commands below is important (see README_fastq-rename.md).

# Check if the 'rename' command is available
if ! command -v rename &> /dev/null; then
    echo "Error: 'rename' command not found. Please install it and try again."
    exit 1
fi


rename 36_R C6_M_NR_R *.gz
rename 35_R C6_M_Red_R *.gz
rename 34_R C5_M_NR_R *.gz
rename 33_R C5_M_Red_R *.gz
rename 32_R A6_F_NR_R *.gz
rename 31_R A6_F_Red_R *.gz
rename 30_R A6_M_NR_R *.gz
rename 29_R A6_M_Red_R *.gz
rename 28_R A5_F_NR_R *.gz
rename 27_R A5_F_Red_R *.gz
rename 26_R A5_M_NR_R *.gz
rename 25_R A5_M_Red_R *.gz
rename 24_R C4_M_NR_R *.gz
rename 23_R C4_M_Red_R *.gz
rename 22_R C3_M_NR_R *.gz
rename 21_R C3_M_Red_R *.gz
rename 20_R A4_F_NR_R *.gz
rename 19_R A4_F_Red_R *.gz
rename 18_R A4_M_NR_R *.gz
rename 17_R A4_M_Red_R *.gz
rename 16_R A3_F_NR_R *.gz
rename 15_R A3_F_Red_R *.gz
rename 14_R A3_M_NR_R *.gz
rename 13_R A3_M_Red_R *.gz
rename 12_R C2_M_NR_R *.gz
rename 11_R C2_M_Red_R *.gz
rename 10_R C1_M_NR_R *.gz
rename 9_R C1_M_Red_R *.gz
rename 8_R A2_F_NR_R *.gz
rename 7_R A2_F_Red_R *.gz
rename 6_R A2_M_NR_R *.gz
rename 5_R A2_M_Red_R *.gz
rename 4_R A1_F_NR_R *.gz
rename 3_R A1_F_Red_R *.gz
rename 2_R A1_M_NR_R *.gz
rename 1_R A1_M_Red_R *.gz

# Indicate completion
echo "Renaming completed successfully."


































