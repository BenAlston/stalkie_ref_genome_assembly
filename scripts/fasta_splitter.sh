#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --mem=6G
#SBATCH -c 6
#SBATCH --time=96:00:00
#SBATCH -e reports/blast/error_fasta.txt
#SBATCH -o reports/blast/output_fasta.txt

#########################################################################################
#       Script Name: fasta.splitter
#       Description: splits a fasta into 500 smaller files
#       Author:      Ben Alston
#       Date:        April 2024
#########################################################################################


module load Anaconda3/2022.05
source activate blast # conda env contains the package "fasta splitter"


species=whitei
assembly=2
wd=/mnt/parscratch/users/bip23bta/ref_genomes/$species
blastdb=/mnt/parscratch/users/bip23bta/ref_genomes/blast_dbs/blast_nt_db/nt
blast='apptainer exec ~/blast_latest.sif'
INPUT_ASSEMBLY=$wd/02-hifiasm/${assembly}_hifiasm_output/*_primary.fa

mkdir $wd/03-QC/blast/${assembly}_primary
cd $wd/03-QC/blast/${assembly}_primary

#
ln -s $INPUT_ASSEMBLY

# splits a given fasta file into 500 smaller ones
fasta-splitter --n-parts 500 $(basename $INPUT_ASSEMBLY) --nopad --out-dir split
