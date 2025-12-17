#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --mem=180G
#SBATCH -c 32
#SBATCH --time=48:00:00
#SBATCH -o reports/output.hifiasm.meigenii_4.txt

#########################################################################################
#       Script Name: Hifiasm
#       Description: genome assembly with hifiasm
#       Author:      Ben Alston
#       Date:        APR 2024 
#########################################################################################

# uses Hifiasm v0.16.1-r375
module load Anaconda3
source activate hifiasm

# Variables
SPECIES=meigenii
ASSEMBLY=4
SEX=Female
WD=/mnt/parscratch/users/bip23bta/ref_genomes/$SPECIES/02-hifiasm/${ASSEMBLY}_omni-c_hifiasm_output_final

# files
HIFI_READS=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/data/long_read/200437_${ASSEMBLY}-Cell?/*.fastq.gz
HI_C_READS_1=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/data/omni_c/*-${SEX}_R1_001.fastq.gz
HI_C_READS_2=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/data/omni_c/*-${SEX}_R2_001.fastq.gz

# ------------ script ------------ #
mkdir -p $WD 
cd $WD 

# symbolic links to input reads - for some reason it wont run if the input reads aren't in the wd
bash -c "ln -fs ${HIFI_READS} ./"
bash -c "ln -fs ${HI_C_READS_1} ./"
bash -c "ln -fs ${HI_C_READS_2} ./"

# hifiasm ----------

hifiasm -o ${SPECIES}_${ASSEMBLY}_omni-c.asm -t 32 --h1 $HI_C_READS_1 --h2 $HI_C_READS_2 $HIFI_READS
