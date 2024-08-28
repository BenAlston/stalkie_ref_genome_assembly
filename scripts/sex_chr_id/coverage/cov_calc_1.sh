#!/bin/bash
#SBATCH --job-name=cov_calc
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=96:00:00
#SBATCH -o reports/cov_calc.output.%J.txt
#SBATCH --array 0-4%5 # change depending on how many tasks you need doing

# load bedtools/2.26.0-0 samtools/1.3.1-0
module load Anaconda3
source activate cov_calc

#########################################################################################
#       Script Name: cov_calc_1.sh
#       Description: First stage in calculating coverage across 5kb genomic windows.
#                    Run as an array with a file of file name prefixes. Takes the .bam file 
#                    following alignment and outputs a sorted, merged, .bed file
#       Author:      Ben Alston
#       Date:        Jun 2024
#########################################################################################
# set wd
cd /mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/males

# read in a list of file prefixes corresponding to the samples I want to map, where each line is a new sample
FILEPREFIXES=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/males/male_fileprefixes.txt

# read this as an array (a specific kind of unix object, kinda like vectors in R)
readarray -t array < $FILEPREFIXES

# Specify that each slurm array job corresponds to a different sample
FILENAME="${array[$SLURM_ARRAY_TASK_ID]}"

echo "Calculating Coverage for $FILENAME"

# sort and index the .bam file
samtools sort ${FILENAME}.bam -o ${FILENAME}_sorted.bam
samtools index ${FILENAME}_sorted.bam

# generates .bed file from .bam
bedtools bamtobed -i ${FILENAME}_sorted.bam > ${FILENAME}_sorted.bam.bed

# merges all features of the .bed file into one
bedtools merge -i ${FILENAME}_sorted.bam.bed > ${FILENAME}_sorted_merged.bed
