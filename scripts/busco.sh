#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --mem=10G
#SBATCH -c 20
#SBATCH --time=24:00:00
#SBATCH -e reports/error.%j_busco.txt
#SBATCH -o reports/output.%j_busco.txt

#########################################################################################
#       Script Name: busco.sh
#       Description: busco assesment for hifiasm assemblies
#       Author:      Ben Alston
#       Date:        Feb 2024
#########################################################################################

# BUSCO/5.6.1
busco ='apptainer exec /users/bip23bta/busco_v5.6.1_cv1.sif busco'

# Paths
export PATH="/path/to/AUGUSTUS/augustus-3.2.3/bin:$PATH"
export PATH="/path/to/AUGUSTUS/augustus-3.2.3/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/path/to/AUGUSTUS/augustus-3.2.3/config/"

# variables
GENOME=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/03-QC/purge_dups/1_primary/round_2/purged.fa
OUT=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/03-QC/purge_dups/1_primary/round_2/

mkdir $OUT

$busco
  -i ${GENOME} \
  -l diptera_odb10 \
  -o /$OUT \
  -m genome \
  -f \
-c 20

