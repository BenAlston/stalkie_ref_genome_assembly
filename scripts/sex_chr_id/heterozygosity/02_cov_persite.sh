#!/bin/bash
#SBATCH --job-name=persite_cov
#SBATCH --mem=20G
#SBATCH -c 2
#SBATCH --time=02:00:00
#SBATCH -o reports/cov_persite.txt

#########################################################################################
#       Script Name: 02_cov_persite.sh
#       Description: calculate coverage per site across the genome for each sample
#                    used in het calculations 
#       Author:      Ben Alston                                                 
#       Date:        Sep 2025                                                    
#########################################################################################

# we need per site coverage for each sample (from less stringent mapped dal_7 bams)
BAMS=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/${1}_bams.txt
# genome file of scaffold names and lengths
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/persite_cov

mkdir -p $WD 
cd $WD 
samtools='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/samtools_1.21--h50ea8bc_0.sif samtools'
$samtools depth -a -f $BAMS -H > ${1}_cov.tsv

#manually edit file prefix because way ive done this is stupid

sed -i 's|/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/01_bams/||g' ${1}_cov.tsv
