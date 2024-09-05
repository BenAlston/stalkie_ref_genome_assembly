#!/bin/bash
#SBATCH --job-name=bowtie_index
#SBATCH --mem=12G
#SBATCH -c 4
#SBATCH --time=12:00:00
#SBATCH -o reports/bowtie_index.txt

#########################################################################################
#       Script Name: bowtie2_index.sh
#       Description: generates index files from a ref assembly for alignment of reads using bowtie2
#       Author:      Ben Alston
#       Date:        May 2024
#########################################################################################

# load bowtie2/2.5.3
module load Anaconda3
source activate bowtie_align

# sample variables
SPECIES=whitei
ASSEMBLY=2

#filepaths
OUTPUT_DIR=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/index_2
REF=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/02-hifiasm/2_hifiasm_output/2_primary.fa

mkdir $OUTPUT_DIR
cd $OUTPUT_DIR

bowtie2-build $REF ${SPECIES}_${ASSEMBLY}
