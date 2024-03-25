#!/bin/bash
#SBATCH --job-name=purge_dups
#SBATCH --mem=16G
#SBATCH -c 12
#SBATCH --time=24:00:00
#SBATCH -e reports/error_purge_dups.txt
#SBATCH -o reports/output_purge_dups.txt

#########################################################################################
#	Script Name:
#	Description:
#	Author:      Ben Alston
#	Date:
#########################################################################################

# load package
module load Anaconda3/2022.05
source activate purge_dups

wd=/mnt/parscratch/users/bip23bta/ref_genomes
output_dir=whitei/03-QC/purge_dups/1_config
cd $wd

~/purge_dups/scripts/pd_config.py \
-l $output_dir/ \
-s whitei/data/long_read/*_1-*/*.fastq.gz \
-n $output_dir/1_config.json \
whitei/02-hifiasm/1_hifiasm_output/1_primary.fa
