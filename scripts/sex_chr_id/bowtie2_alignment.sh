#!/bin/bash
#BATCH --job-name=bowtie2_align
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=96:00:00
#SBATCH -e reports/bowtie_align.error.%j.txt
#SBATCH -o reports/bowtie_align.output.%j.txt

#########################################################################################
#       Script Name: bowtie2_align.sh
#       Description: Aligns short reads to female ref, outputing a .sam file
#       Author:      Ben Alston
#       Date:        May 2024
#########################################################################################

# load bowtie2/2.4.1 and samtools/1.13
module load Anaconda3
source activate bowtie_align

WD=/mnt/parscratch/users/bip23bta/ref_genomes/
READS_1=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Trimmed/female/*R1_001.fastq.gz
READS_2=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Trimmed/female/*R2_001.fastq.gz
REF_INDEX=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/index_1/whitei_1

cd $WD 

bowtie2 --mp 10000 -x $REF_INDEX -1 $READS_1 -2 $READS_2 | samtools view -bS - > whitei/04-sex_chromosomes/female_whitei.bam 
# high mismatch penalty (mp) effectivley removes mismatches
