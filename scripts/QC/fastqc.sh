#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH -e reports/fastqc_error.txt
#SBATCH -o reports/fastqc_output.txt

#########################################################################################
#       Script Name: fastqc
#       Description: FastQC Quality Assessment for Raw Data
#       Author:      Ben Alston
#       Date:        Jan 2024
#########################################################################################

module load Anaconda3
source activate qc

# fastqc/0.12.1 multiqc/1.0


READS=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Trimmed_trimmomatic/*.paired.fastq
OUTPUT=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Trimmed_trimmomatic/
 
cd /mnt/parscratch/users/bip23bta/ref_genomes 

for FILE in $READS    
 
do fastqc $FILE -o $OUTPUT 

done
