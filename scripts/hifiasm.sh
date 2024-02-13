#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --mem=36G
#SBATCH -c 32
#SBATCH --time=96:00:00
#SBATCH -e error%j_hifiasm.txt
#SBATCH -o output.%j_hifiasm.txt
#SBATCH --array 0-1%2 # change depending on how many assemblies you're making, needs to start with 0

#########################################################################################
#	Script Name: hifiasm
#	Description: genome assembly with hifiasm, seperate for M & F in each species (plus a whitei driver F)
#	Author:      Ben Alston
#	Date:        Feb 2024
#########################################################################################

# load package
module load Anaconda3
source activate hifiasm

# variables
wd=/mnt/parscratch/users/bip23bta/ref_genomes
species=dalmanni
samples=(6 7)
assembly=${samples[$SLURM_ARRAY_TASK_ID]} # make the array task ids correspond with the samples
data=${wd}/$species/data/long_read/200437_${assembly}-Cell?/*.fastq.gz # input data
output_log=${wd}/$species/02-hifiasm/${SLURM_ARRAY_TASK_ID}_${assembly}.output.txt
# these correspond to the liverpool prefixes, matched to our ones in <shared lab folder>/2024_Stalkie_PacBio_HiFI/>

# ------------ script
cd ${wd}

# make and set dirs
echo starting $assembly >> $output_log

echo mkdir ${wd}/$species/02-hifiasm/${assembly}_hifiasm_output >> ${output_log}
mkdir ${wd}/$species/02-hifiasm/${assembly}_hifiasm_output

echo cd ${wd}/$species/02-hifiasm/${assembly}_hifiasm_output >> ${output_log}
cd ${wd}/$species/02-hifiasm/${assembly}_hifiasm_output

# symbolic link(s)
echo bash -c "ln -s ${data} ./" >> ${output_log}
bash -c "ln -s ${data} ./"

# hifiasm
echo hifiasm -o ${species}_${assembly}.asm -t 32 *.fastq.gz >> ${output_log}
hifiasm -o ${species}_${assembly}.asm -t 32 *.fastq.gz

# Give your script more cores
# It's not happy with just one
# It's doing it's best
