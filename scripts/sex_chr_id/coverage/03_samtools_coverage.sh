#!/bin/bash
#SBATCH --job-name=cov_calc
#SBATCH --mem=4G
#SBATCH -c 2
#SBATCH --time=24:00:00
#SBATCH -o reports/cov_calc/%J.cov.output.txt
#SBATCH --array 0-4%5 

# load bedtools/2.26.0-0 samtools/1.3.1-0
module load Anaconda3
source activate cov_calc

#########################################################################################
#       Script Name: 02_cov_calc_1.sh
#       Description: calculates coverage stats per sample per contig
#       Author:      Ben Alston
#       Date:        Nov 2025
#########################################################################################

# array ----
# read in a list of file prefixes corresponding to the samples I want to map, where each line is a new sample
FILEPREFIXES=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/${1}_sample_list.txt
# read this as an array (a specific kind of unix object, kinda like vectors in R)
readarray -t array < $FILEPREFIXES
# Specify that each slurm array job corresponds to a different sample
FILENAME="${array[$SLURM_ARRAY_TASK_ID]}"


samtools='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/samtools_1.21--h50ea8bc_0.sif samtools'


# filepaths ----
REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/refs/dal_7_scaffolded.fa
BAM=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/bams/${FILENAME}_sorted.dedup.bam
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/cov_per_contig/${1}
echo "Calculating Coverage for $FILENAME"

mkdir -p $WD 
cd $WD 

$samtools faidx $REF
cut -f1,2 ${REF}.fai | awk '{print $1":1-"$2}' > chrsizes.list
cat chrsizes.list | while read region
do
        echo $region
        $samtools coverage --no-header -r $region ${BAM} >> ${FILENAME}_cov.tsv
done
