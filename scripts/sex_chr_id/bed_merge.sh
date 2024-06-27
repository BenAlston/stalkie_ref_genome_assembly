#!/bin/bash
#SBATCH --job-name=intersect_bed
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=96:00:00
#SBATCH -o reports/bed.merge.txt

#########################################################################################
#       Script Name: bed_merge.sh
#       Description: merges all individual .bed files into one then seperates this into 5kb genomic windows.                                
#       Author:      Ben Alston                                                 
#       Date:        Jun 2024                                                    
#########################################################################################

# load bedtools/2.26.0-0 samtools/1.3.1-0
module load Anaconda3
source activate cov_calc

# set wd
cd /mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/

bedtools multiinter -i females/*_sorted.bam.bed males/*_sorted.bam.bed > males/merged_all.sorted.bed

bedtools merge -i males/merged_all.sorted.bed > males/merged_all.sorted.2.bed

bedtools makewindows -b males/merged_all.sorted.2.bed -w 5000 > males/merged_all.windows.2.bed

#cp males/merged_all.windows.bed  females/merged_all.windows.bed
