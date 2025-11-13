#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=96:00:00
#SBATCH -e reports/trimmomatic.error.%J.txt
#SBATCH -o reports/trimmomatic.output.%J.txt

#########################################################################################
#       Script Name: trim.sh
#       Description: Trim reads with trimmomatic                                
#       Author:      Ben Alston                                                 
#       Date:        Jun 2024                                                    
#########################################################################################

module load Anaconda3
source activate trimmomatic

# filepaths
DIR=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Raw/*
OUTPUT=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Trimmed_trimmomatic
cd /mnt/parscratch/users/bip23bta/ref_genomes

for READS in $DIR
do
trimmomatic PE \
$READS/*R1_001.fastq.gz $READS/*R2_001.fastq.gz \
${OUTPUT}/$(basename $READS/*R1_001.fastq.gz).paired.fastq ${OUTPUT}/$(basename $READS/*R1_001.fastq.gz).unpaired.fastq \
${OUTPUT}/$(basename $READS/*R2_001.fastq.gz).paired.fastq ${OUTPUT}/$(basename $READS/*R2_001.fastq.gz).unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
done
