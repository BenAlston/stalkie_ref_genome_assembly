#!/bin/bash
#SBATCH --job-name=jellyfish
#SBATCH --mem=24G
#SBATCH --time=96:00:00
#SBATCH -e jellyfish_error.txt
#SBATCH -o jellyfish_output.txt

#########################################################################################
#	Script Name: jellyfish
#	Description: kmer counting with jellyfish
#	Author:      Ben Alston
#	Date:        Feb 2024
#########################################################################################

# jellyfish/2.2.10
module load Anaconda3
source activate jellyfish

# variables
wd=/mnt/parscratch/users/bip23bta/ref_genomes/
output=01-QC/long_read/jellyfish
species=whitei
samples="1 2 3"

# loop
for assembly in $samples

do

cd ${wd}/$species
mkdir $output/${assembly}_output

jellyfish count <(zcat data/long_read/200437_${assembly}-Cell?/*.fastq.gz) -C -m 21 -s 1G -t 25 -o ${output}/${assembly}_output/$species.jf
jellyfish histo ${output}/${assembly}_output/$species.jf > ${output}/${assembly}_output/$species_${assembly}.histo

done
