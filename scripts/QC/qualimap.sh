#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --mem=100G
#SBATCH -c 16
#SBATCH --time=96:00:00
#SBATCH -o reports/output.%j_qualimap.txt

#########################################################################################
#       Script Name: qualimap
#       Descrption:  aligment with minimap2 and quality assesment with qualimap
#       Author:      Ben Alston
#       Date:        Jul 2024
#########################################################################################

# load minimap2/2.17-r941, qualimap/v.2.2.2
module load Anaconda3/2022.05
source activate qualimap 

# input variables
SPECIES=dalmanni
ASSEMBLY=6

#filepaths
WD=/mnt/parscratch/users/bip23bta/ref_genomes/$SPECIES
REF_GENOME=02-hifiasm/${ASSEMBLY}_hifiasm_output/${ASSEMBLY}_primary.fa
READS=data/long_read/*_${ASSEMBLY}-*/*.fastq.gz # dir containing long read data to be remapped
OUTPUT=02-hifiasm/${ASSEMBLY}_remapping

# ----- script
cd $WD

echo running for $SPECIES $ASSEMBLY 

# map to ref with minimap
minimap2 -t 16 -a $REF_GENOME $READS > $OUTPUT/${ASSEMBLY}_primary.bam

# qualimap requires a sorted bam file
samtools sort $OUTPUT/${ASSEMBLY}_primary.bam -o $OUTPUT/${ASSEMBLY}_primary_sorted.bam

# remove uneeded unsorted bam file
rm $OUTPUT/${ASSEMBLY}_primary.bam

# run qualimap, this step requires lots of memory
qualimap bamqc -bam $OUTPUT/${ASSEMBLY}_primary_sorted.bam --java-mem-size=100G -outdir $OUTPUT/qualimap_results/
