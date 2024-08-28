#!/bin/bash
#SBATCH --job-name=cov_calc2
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=96:00:00
#SBATCH -o reports/cov_calc2.output.%J.txt
#SBATCH --array 0-4%5 # change depending on how many tasks you need doing

#########################################################################################
#       Script Name: cov_calc_2.sh
#       Description: 2nd stage of calculating coverage for 5kb genomic windows. Needs a merged, sorted 
#                    .bed file of all files being compared, that has been split into 5kb windows using bed_merged.sh.  
#                    Outputs a .tsv containing 4 cols: contig_name, window_start, window_end, coverage (in that order)
#       Author:      Ben Alston                                                 
#       Date:        Jun 2024                                                    
#########################################################################################

# load bedtools/2.26.0-0 samtools/1.3.1-0
module load Anaconda3
source activate cov_calc


# set wd
cd /mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/males

# read in a list of file prefixes corresponding to the samples I want to map, where each line is a new sample
FILEPREFIXES=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/males/male_fileprefixes.txt

# read this as an array (a specific kind of unix object, kinda like vectors in R)
readarray -t array < $FILEPREFIXES

# Specify that each slurm array job corresponds to a different sample
FILENAME="${array[$SLURM_ARRAY_TASK_ID]}"
MERGED_BED=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/males/merged_all.sorted.bed # all should be mapped to the same bed file, otherwise it's impossible to merge male and male reads

echo "Calculating coverage for $FILENAME (part 2)"

# counts the number of alignments that overlap with the intervals specified by the .bed file
bedtools multicov -bams ${FILENAME}_sorted.bam -bed merged_all.windows.bed > ${FILENAME}_m_cov.tsv

