#!/bin/bash
#SBATCH --job-name=quast
#SBATCH --mem=40G
#SBATCH -c 32
#SBATCH --time=24:00:00
#SBATCH -e reports/error%j_quast.txt
#SBATCH -o reports/output.%j_quast.txt
#SBATCH --array 0%1 # change depending on how many assemblies you're making, needs to start with 0


#########################################################################################
#	Script Name: quast
#	Description: basic ref genome summary statistics with quast
#	Author:      Ben Alston
#	Date:        Feb 2024
#########################################################################################

module load Anaconda3/2022.05 # specific version of anaconda needed on the hpc atm
source activate quast

# variables
wd=/mnt/parscratch/users/bip23bta/ref_genomes

species=dalmanni
assembly=7 # make the array task ids correspond with the samples
data=02-hifiasm/${assembly}_hifiasm_output # input data
#output_log=02-hifiasm/jellyfish/${SLURM_ARRAY_TASK_ID}_${assembly}.output.txt 

cd $wd/$species

equast $data/${species}_${assembly}.asm.bp.p_ctg.fa -o $data/quast_output -t 32





