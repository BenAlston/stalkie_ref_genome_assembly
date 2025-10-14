#!/bin/bash
#BATCH --job-name=genomicsDB
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH --time=96:00:00
#SBATCH -o reports/03_%J_genomicsDBimport.txt

#########################################################################################
#       Script Name: 03_genomicsDBimport.sh
#       Description: combines vcf files from a list of sample bam files
#                    input is a .map file of vcf direct filepaths
#       Author:      Ben Alston
#       Date:        May 2024
#########################################################################################

# load gatk/v4.3.0.0
module load GATK

# working dir
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/03_genomicsdb

# filenames and output
REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/refs/dal_6_blobtools_scaffolded.fa
OUTPUT=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/03_genomicsdb/${1}_db
MAP=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/${1}_samples.map
INTERVALS=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/dal_6_intervals.list

# ---- WD ---- #
mkdir -p $WD 
cd $WD 

# ---- run genomicsdbimport ---- #
gatk GenomicsDBImport \
        --genomicsdb-workspace-path $OUTPUT \
        --sample-name-map $MAP \
        --intervals $INTERVALS
# --intervals is required even when running on whole genome. Here dal_7_contigs.list is just a list of all dal 7 contig names
