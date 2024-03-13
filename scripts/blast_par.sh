#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --mem=10G
#SBATCH -c 16
#SBATCH --time=96:00:00
#SBATCH -e reports/error_par_blast.txt
#SBATCH -o reports/output_par_blast.txt
#SBATCH --array 1-500%50 # concurrent job limit is 50, othewise i'd do 500

#########################################################################################
#      Script Name: blast_par
#      Description: runs blast on a genome assembly subdivided into 500 smaller .fa files
#      Author:      Ben Alston
#      Date: March 2024
#########################################################################################

wd=/mnt/parscratch/users/bip23bta/ref_genomes
species=dalmanni
assembly=6
input_file=${assembly}_primary.fa

cd $wd/$species/03-QC/blast

# using the apptainer image of blast since the other versions dont work
apptainer exec ~/blast_latest.sif \
        blastn -db ${wd}/blast_nt_db/nt \
       -query split/${assembly}_primary.part-${SLURM_ARRAY_TASK_ID}.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -out out/${SLURM_ARRAY_TASK_ID}_blast.out
