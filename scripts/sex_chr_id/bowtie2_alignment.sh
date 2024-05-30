#!/bin/bash
#SBATCH --job-name=bowtie2_alignment
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH -e reports/bowtie_align.%J.txt
#SBATCH -o reports/bowtie_align.%J.txt

#########################################################################################
#	Script Name: bowtie2_alignment.sh
#	Description: Aligns short reads to female ref
#	Author:      Ben Alston
#	Date:        May 2024
#########################################################################################

# load bowtie2/2.5.3
bowtie2='apptainer exec /users/bip23bta/bowtie2_2.5.3--6dd1a06.sif bowtie2'

# sample variables
SPECIES=whitei
ASSEMBLY=1

#filepaths
WD=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/
# Selects only female reads, forward and reverse are seperate (1 & 2)
READS_1=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Trimmed/Sample_??-{A01,A02,A03,A04,A05}/*R1_001.fastq.gz
READS_2=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/data/short_read/Trimmed/Sample_??-{A01,A02,A03,A04,A05}/*R2_001.fastq.gz
# select reference from purge dups output
REF_INDEX=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/index_1/whitei_1

cd $WD

bowtie2 --mp 10000 -x $REF_INDEX -1 $READS_1 -2 $READS_2 -S 04-sex_chromosomes/whitei_1_f.sam
# high mismatch penalty (mp) effectivley removes mismatches





