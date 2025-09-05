#!/bin/bash
#SBATCH --job-name=bed_merge
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=24:00:00
#SBATCH -o reports/bed.merge.txt

#########################################################################################
#       Script Name: 03_bed_merge.sh
#       Description: merges all individual .bed files into one then seperates this into 5kb genomic windows.                                
#       Author:      Ben Alston                                                 
#       Date:        Jun 2024                                                    
#########################################################################################

# load bedtools/2.26.0-0 samtools/1.3.1-0
module load Anaconda3
source activate cov_calc

WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/

# set wd
cd $WD 

bedtools multiinter -i F/*.bed M/*.bed > all.bed

bedtools merge -i all.bed > all.merged.bed

bedtools makewindows -b all.merged.bed -w 5000 > all.merged.5kbwindows.bed
