#!/bin/bash
#SBATCH --job-name=hifiasm_whitei
#SBATCH --mem=180G
#SBATCH -c 48
#SBATCH --time=96:00:00
#SBATCH -e reports/error%j_hifiasm.txt
#SBATCH -o reports/output.%j_hifiasm.txt


#########################################################################################
#	Script Name: hifiasm
#	Description: genome assembly with hifiasm
#	Author:      Ben Alston
#	Date:        April 2024
#########################################################################################

# uses Hifiasm v0.16.1-r375
module load Anaconda3
source activate hifiasm

# Variables
WD=/mnt/parscratch/users/bip23bta/ref_genomes
SPECIES=dalmanni
ASSEMBLY=7
DATA=${WD}/$SPECIES/data/long_read/*${ASSEMBLY}-Cell?/*.fastq.gz

# ------------ script ------------ #
cd ${WD}/$SPECIES/02-hifiasm

# make and set dirs -----------
mkdir ${WD}/$SPECIES/02-hifiasm/${ASSEMBLY}_hifiasm_output

cd ${WD}/$SPECIES/02-hifiasm/${ASSEMBLY}_hifiasm_output

# symbolic links to input reads - for some reason it wont run if the input reads aren't in the wd
bash -c "ln -s ${DATA} ./"

# hifiasm ----------
hifiasm -o ${SPECIES}_${ASSEMBLY}.asm -t 48 $DATA

# convert .gfa to .fa
awk '/^S/{print ">"$2"\n"$3}' ${SPECIES}_${ASSEMBLY}.asm.bp.p_ctg.gfa  | fold > ${ASSEMBLY}_primary.fa
