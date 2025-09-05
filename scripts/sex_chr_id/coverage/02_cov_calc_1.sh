#!/bin/bash
#SBATCH --job-name=cov_calc
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=24:00:00
#SBATCH -o reports/cov_calc.output.txt
#SBATCH --array 0-4%5 

# load bedtools/2.26.0-0 samtools/1.3.1-0
module load Anaconda3
source activate cov_calc

#########################################################################################
#       Script Name: 02_cov_calc_1.sh
#       Description: First stage in calculating coverage across 5kb genomic windows.
#                    Run as an array with a file of file name prefixes. Takes the .bam file 
#                    following alignment and outputs a sorted, merged, .bed file
#       Author:      Ben Alston
#       Date:        Jun 2024
#########################################################################################

# array ----
# read in a list of file prefixes corresponding to the samples I want to map, where each line is a new sample
FILEPREFIXES=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/${1}_sample_list.txt
# read this as an array (a specific kind of unix object, kinda like vectors in R)
readarray -t array < $FILEPREFIXES
# Specify that each slurm array job corresponds to a different sample
FILENAME="${array[$SLURM_ARRAY_TASK_ID]}"

# filepaths ----
BAM=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/${FILENAME}_sorted.dedup.bam
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/${1}

echo "Calculating Coverage for $FILENAME"

mkdir -p $WD 
cd $WD 

# generates .bed file from .bam
bedtools bamtobed -i ${BAM} > ${FILENAME}.bed

# merges all features of the .bed file into one
bedtools merge -i ${FILENAME}.bed > ${FILENAME}.merged.bed
