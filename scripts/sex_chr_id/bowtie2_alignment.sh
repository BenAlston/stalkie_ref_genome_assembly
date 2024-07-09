#!/bin/bash
#BATCH --job-name=bowtie2_align
#SBATCH --mem=16G
#SBATCH -c 32
#SBATCH --time=96:00:00
#SBATCH -o reports/bowtie_align.output.%j.txt
#SBATCH --array 0-4%5 # change depending on how many tasks you need doing

#########################################################################################
#       Script Name: bowtie2_align.sh
#       Description: Aligns illumina paired reads to a ref. Swap between males and females with find and replace
#       Author:      Ben Alston
#       Date:        May 2024
#########################################################################################

# load bowtie2/2.5.3 and samtools
module load Anaconda3
source activate bowtie_align

WD=/mnt/parscratch/users/bip23bta/ref_genomes/
FILEPREFIXES=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/females/female_fileprefixes.txt

readarray -t array < $FILEPREFIXES # a file of line seperated file prefixes

FILENAME="${array[$SLURM_ARRAY_TASK_ID]}"

READS_1=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Trimmed/female/??-${FILENAME}*R1_001.fastq.gz
READS_2=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Trimmed/female/??-${FILENAME}*R2_001.fastq.gz
REF_INDEX=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/index_1/whitei_1

echo "running bowtie2 on sample $FILENAME"

cd $WD 

bowtie2 --threads 32 --mp 10000 -x $REF_INDEX -1 $READS_1 -2 $READS_2 | samtools view -bS - > whitei/04-sex_chromosomes/females/${FILENAME}.bam
# high mismatch penalty (mp) effectivley removes mismatches
