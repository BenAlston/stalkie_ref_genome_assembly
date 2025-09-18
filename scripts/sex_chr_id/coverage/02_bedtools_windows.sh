#!/bin/bash
#SBATCH --job-name=intersect_bed
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH --time=02:00:00
#SBATCH -o reports/bed_windows.txt

#########################################################################################
#       Script Name: bed_windows.sh
#       Description: make bed file of n kb windows, required by bedtools multicov, to specify windows to generate coverage across
#       Author:      Ben Alston                                                 
#       Date:        Sep 2025                                                    
#########################################################################################

# load bedtools/2.26.0-0
module load Anaconda3
source activate cov_calc

# --- filepaths ---- #
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/
REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/refs/dal_7_scaffolded.fa
GENOME=dal_7_chrsizes.g

# set wd
cd $WD 

# make a GENOME file containing scaffold names and sizes from reference fasta (a required input of bedtools makewindows)
cat ../refs/dal_7_scaffolded.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $GENOME

# use this file to generate a bed file of 100kb windows in which to calculate coverage
bedtools makewindows -g $GENOME  -w 100000 > all.merged.100kbwindows.bed
