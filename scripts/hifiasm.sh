#!/bin/bash
#SBATCH --job-name=hifiasm_whitei
#SBATCH --mem=180G
#SBATCH -c 48
#SBATCH --time=96:00:00
#SBATCH -e reports/error%j_hifiasm.txt
#SBATCH -o reports/output.%j_hifiasm.txt
#SBATCH --array 0%2 # change depending on how many assemblies you're making, needs to start with 0

#########################################################################################
#	Script Name: hifiasm
#	Description: genome assembly with hifiasm, seperate for M & F in each species (plus a whitei driver F)
#	Author:      Ben Alston
#	Date:        Feb 2024
#########################################################################################

# load package
module load Anaconda3
source activate hifiasm

# Variables
wd=/mnt/parscratch/users/bip23bta/ref_genomes
species=whitei
samples=(1 2 3)
assembly=${samples[$SLURM_ARRAY_TASK_ID]}
data=${wd}/$species/data/long_read/*${assembly}-Cell?/*.fastq.gz
output_log=${wd}/$species/02-hifiasm/${SLURM_ARRAY_TASK_ID}_${assembly}.output.txt
# these correspond to the liverpool prefixes, matched to our ones in <shared lab folder>/2024_Stalkie_PacBio_HiFI/>

# ------------ script
cd ${wd}

# make and set dirs
date
echo making dirs and symlinks
mkdir ${wd}/$species/02-hifiasm/${assembly}_hifiasm_output

cd ${wd}/$species/02-hifiasm/${assembly}_hifiasm_output

# symbolic links

bash -c "ln -s ${data} ./"

# hifiasm
date
echo starting hifiasm
hifiasm -o ${species}_${assembly}.asm -t 48 *.fastq.gz
date

# convert gfa to fa
awk '/^S/{print ">"$2"\n"$3}' ${species}_${assenbly}.asm.bp.p_ctg.gfa  | fold > ${assembly}_primary.fa
