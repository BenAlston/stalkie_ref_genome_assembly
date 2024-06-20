#!/bin/bash
#SBATCH --job-name=cov_calc
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=96:00:00
#SBATCH -e reports/cov_calc.error.%J.txt
#SBATCH -o reports/cov_calc.output.%J.txt

#########################################################################################
#       Script Name: cov_calc.sh
#       Description: cov                                
#       Author:      Ben Alston                                                 
#       Date:        Jun 2024                                                    
#########################################################################################

# load bedtools/2.26.0-0 samtools/1.3.1-0
module load Anaconda3
source activate cov_calc

cd /mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/

# sort and index the .bam file
samtools sort female_whitei.bam -o female_whitei_sorted.bam
samtools index female_whitei_sorted.bam

# merges all features of the .bed file into one
# will need to check if this is the right way of doing it
bedtools merge -i female_whitei_sorted.bam.bed > female_whitei_sorted_merged.bed

#bedtools bamtobed -i female_whitei_sorted.bam | bedtools makewindows -b - -w 5000 -s 5000 > female_whitei_windows.bed

# generates .bed file from .bam
bedtools bamtobed -i female_whitei_sorted.bam > female_whitei_sorted.bam.bed

bedtools makewindows -b female_whitei_sorted_merged.bed -w 5000

# counts the number of alignments that overlap with the intervals specified by the .bed file
bedtools multicov -bams female_whitei_sorted.bam -bed female_whitei_windows.bed > female_cov.csv



