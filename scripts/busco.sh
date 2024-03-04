#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --mem=10
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

# had to use a docker image instead of a conda env to get most recent version, this needs to be called before the comman

export PATH="/path/to/AUGUSTUS/augustus-3.2.3/bin:$PATH"
export PATH="/path/to/AUGUSTUS/augustus-3.2.3/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/path/to/AUGUSTUS/augustus-3.2.3/config/"



# variables
wd=/mnt/parscratch/users/bip23bta/ref_genomes

species=dalmanni
genome=7 # make the array task ids correspond with the samples
data_path=02-hifiasm/${genome}_hifiasm_output # input data
output=${data_path}
output_log=02-hifiasm/jellyfish/${SLURM_ARRAY_TASK_ID}_${genome}.output.txt
database=diptera_odb10

cd $wd/$species

# call docker image before busco command
apptainer exec /users/bip23bta/busco_v5.6.1_cv1.sif \
busco -i ${data_path}/${species}_${genome}.asm.bp.p_ctg.fa \
-l $database \
-o ${output}/BUSCO_out_${genome} \
-m genome \
-f \
-c 20
# -f=fource, overrides previous output dir
# -c=cores
