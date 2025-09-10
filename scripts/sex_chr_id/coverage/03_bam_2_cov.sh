#!/bin/bash
#SBATCH --job-name=bam_2_cov
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=12:00:00
#SBATCH -o reports/bam_2_cov_output.txt
#SBATCH --array 0-4%5

#########################################################################################
#       Script Name: 04_bam_2_cov.sh
#       Description: Extracts coverage info from a bam file using intervals specified by a bed file
#       Author:      Ben Alston                                                 
#       Date:        Jun 2024                                                    
#########################################################################################

# load bedtools/2.26.0-0 samtools/1.3.1-0
module load Anaconda3
source activate cov_calc

# array ----
# read in a list of file prefixes corresponding to the samples I want to map, where each line is a new sample
FILEPREFIXES=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/${1}_sample_list.txt
# read this as an array (a specific kind of unix object, kinda like vectors in R)
readarray -t array < $FILEPREFIXES
# Specify that each slurm array job corresponds to a different sample
FILENAME="${array[$SLURM_ARRAY_TASK_ID]}"

# filepaths ----
BED=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/all.merged.5kbwindows.bed
BAM=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/${FILENAME}_sorted.dedup.bam
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/cov_output

mkdir -p $WD 
cd $WD 

# counts the number of alignments that overlap with the intervals specified by the .bed file
bedtools multicov -bams $BAM -bed $BED > ${FILENAME}_${1}_cov.tsv
