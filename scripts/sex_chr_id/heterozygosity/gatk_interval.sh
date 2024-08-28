#!/bin/bash
#SBATCH --job-name=interval_gatk
#SBATCH --mem=10G
#SBATCH -c 6
#SBATCH --time=24:00:00
#SBATCH -o reports/output.gatk.txt

#########################################################################################
#       Script Name: qualimap
#       Descrption:  
#       Author:      Ben Alston
#       Date:        jul 2024
#########################################################################################

module load GATK # load gatk/4.3.0.0 already installed onto the hpc

#filepath variables
WD=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes
REF=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/02-hifiasm/1_hifiasm_output/1_primary.fa

cd ${WD}/het

gatk SplitIntervals -R ${REF} --scatter-count 60 -O ./
