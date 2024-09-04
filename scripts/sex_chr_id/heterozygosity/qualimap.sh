#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --mem=100G
#SBATCH -c 16
#SBATCH --time=96:00:00
#SBATCH -o reports/output_qualimap.txt

#########################################################################################
#       Script Name: qualimap
#       Descrption:  quality assesment with qualimap
#       Author:      Ben Alston
#       Date:        Jul 2024
#########################################################################################

# load minimap2/2.17-r941, qualimap/v.2.2.2
module load Anaconda3/2022.05
source activate qualimap 

#filepaths
FILES=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/het/female/*.bam
OUTPUT=/mnt/parscratch/users/bip23bta/ref_genomes/whitei/04-sex_chromosomes/het/female

for BAM in $FILES
do
mkdir $OUTPUT/$(basename -s .bam $BAM)_qualimap
# run qualimap
qualimap bamqc -bam $BAM --java-mem-size=100G -outdir $OUTPUT/$(basename -s .bam $BAM)_qualimap
done
