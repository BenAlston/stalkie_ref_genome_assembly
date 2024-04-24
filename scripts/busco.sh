#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --mem=10G
#SBATCH -c 20
#SBATCH --time=24:00:00
#SBATCH -e reports/error.%j_busco.txt
#SBATCH -o reports/output.%j_busco.txt
#SBATCH --array 0%1 # change depending on how many assemblies you're making, needs to start with 0

#########################################################################################
#	Script Name: busco.sh
#	Description: busco assesment for hifiasm assemblies
#	Author:      Ben Alston
#	Date:        Feb 2024
#########################################################################################

# Paths requird by BUSCO, I do not know why, but I shall not question the machine
export PATH="/path/to/AUGUSTUS/augustus-3.2.3/bin:$PATH"
export PATH="/path/to/AUGUSTUS/augustus-3.2.3/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/path/to/AUGUSTUS/augustus-3.2.3/config/"

# variables
wd=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/03-QC/blobtools
species=whitei
GENOME=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/03-QC/purge_dups/1_primary/purged.fa


cd $wd/$species/03-QC/purge_dups/1_primary

apptainer exec /users/bip23bta/busco_v5.6.1_cv1.sif busco -i ${GENOME} \
-l diptera_odb10 \
-o busco_out \
-m genome \
-f \
-c 20

