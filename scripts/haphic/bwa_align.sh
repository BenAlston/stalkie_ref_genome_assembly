#!/bin/bash
#BATCH --job-name=bwa_align
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH --time=96:00:00
#SBATCH -o reports/bwa_align.txt

#########################################################################################
#       Script Name: bwa_align.sh
#       Description: Aligns short reads to ref
#       Author:      Ben Alston
#       Date:        May 2024
#########################################################################################

# load samtools and samblaster
module load Anaconda3
source activate haphic

# load bwa version 0.7.17-r1188
bwa='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/bwa_v0.7.17_cv1.sif bwa'

REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/dalmanni_6/6_omni-c_phasedhaps.fa
READS1=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/data/omni_c/Tdalmanni-Male_R1_001.fastq.gz
READS2=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/data/omni_c/Tdalmanni-Male_R2_001.fastq.gz
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/dalmanni_6

mkdir -p $WD 
cd $WD 

#$bwa index $REF
$bwa mem -t 8 -5SP $REF $READS1 $READS2 | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o tdal_hic_aligned.bam
