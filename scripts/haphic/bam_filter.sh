#!/bin/bash
#BATCH --job-name=bam_filter
#SBATCH --mem=16G
#SBATCH -c 10
#SBATCH --time=48:00:00
#SBATCH -o reports/bam_filter.txt

#########################################################################################
#       Script Name: bam_filter.sh
#       Description: runs haphic bam filter script
#       Author:      Ben Alston
#       Date:        April 2025
#########################################################################################

# load samtools and samblaster
module load Anaconda3
source activate haphic

# filter alignment
BAM=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/dalmanni_6/tdal_hic_aligned.bam
FILTER_BAM=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/HapHiC/utils/filter_bam

cd /mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/dalmanni_6/

$FILTER_BAM $BAM 1 --nm 3 --threads 10 | samtools view - -b -@ 10 -o tdal_hic_filtered.bam
