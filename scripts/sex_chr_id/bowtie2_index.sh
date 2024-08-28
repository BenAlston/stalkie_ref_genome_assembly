#!/bin/bash
#SBATCH --job-name=bowtie_index
#SBATCH --mem=12G
#SBATCH -c 4
#SBATCH --time=24:00:00
#SBATCH -e reports/bowtie_index.txt
#SBATCH -o reports/bowtie_index.txt

#########################################################################################
#	Script Name: bowtie2_index.sh
#	Description: generates index files from a ref assembly for alignment of reads using bowtie2
#	Author:      Ben Alston
#	Date:        May 2024
#########################################################################################

# load bowtie2/2.5.3
bowtie2-build='apptainer exec /users/bip23bta/bowtie2_2.5.3--6dd1a06.sif bowtie2-build'

# sample variables
SPECIES=whitei
ASSEMBLY=1

#filepaths
OUTPUT_DIR=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/index_1
REF=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/02-hifiasm/1_hifiasm_output/1_primary.fa

mkdir $OUTPUT_DIR
cd $OUTPUT_DIR

$bowtie2-build $REF ${SPECIES}_${ASSEMBLY}


