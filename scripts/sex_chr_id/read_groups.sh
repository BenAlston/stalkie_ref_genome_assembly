#!/bin/bash
#BATCH --job-name=readgroups
#SBATCH --mem=40G
#SBATCH -c 16
#SBATCH --time=96:00:00
#SBATCH -o reports/read_groups_output.txt
#SBATCH --array 0-4%5 # change depending on how many tasks you need doing

#########################################################################################
#       Script Name: bowtie2_align.sh
#       Description: Aligns short reads to a ref, outputting a sorted .bam file
#       Author:      Ben Alston
#       Date:        May 2024
#########################################################################################

# load bowtie2/2.5.3 and samtools
module load Anaconda3
source activate bowtie_align

# sample info
SEX=male
ASSEMBLY=1
SPECIES=whitei
REF_INDEX=index_${ASSEMBLY}/${SPECIES}_${ASSEMBLY}
SEXCHR=x

# working dir
WD=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/04-sex_chromosomes

# arrays 
FILEPREFIXES=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/04-sex_chromosomes/${SEX}_fileprefixes.txt
readarray -t array < $FILEPREFIXES
FILENAME="${array[$SLURM_ARRAY_TASK_ID]}"

# reads and index
READS_1=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/data/short_read/Trimmed/${SEX}/??-${FILENAME}*R1_001.fastq.gz
READS_2=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/data/short_read/Trimmed/${SEX}/??-${FILENAME}*R2_001.fastq.gz
REF_INDEX=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/04-sex_chromosomes/index_${ASSEMBLY}/${SPECIES}_${ASSEMBLY}

#############################
# 1 - Alignment     	    #	
#############################

echo "running  on sample $FILENAME"

cd $WD

#bowtie2 --threads 32 -x $REF_INDEX -1 $READS_1 -2 $READS_2 | samtools view -bS | samtools sort -@ 16 -o het/${SEX}/${FILENAME}.sorted.bam



#############################
# 2 - Picard ADD read group #
#############################
module load GATK

mkdir -p het/${SEX}/sorted_dedup_bam/  # Create working directory for merged bam file without duplicates


RGID=$(zcat ${READS_1} | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4) # Extract RGID from read header
RGLB="Solexa"
RGSM="$FILENAME"
RGPU=$RGID.$RGLB
RGPL="Illumina"

gatk AddOrReplaceReadGroups \
 I=het/${SEX}/${FILENAME}.sorted.bam \
 O=het/${SEX}/${FILENAME}.sorted.RG.bam \
 RGID=$RGID \
 RGLB=$RGLB \
 RGPL=$RGPL \
 RGPU=$RGPU \
 RGSM=$RGSM \
 SORT_ORDER=coordinate \
 CREATE_INDEX=TRUE


#############################
# 3 - Picard MarkDuplicates #
#############################
cd $WD
mkdir TMP
gatk MarkDuplicates \
 --REMOVE_DUPLICATES TRUE \
 --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
 --VALIDATION_STRINGENCY SILENT \
 --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
 --TMP_DIR TMP \
 -I het/${SEX}/${FILENAME}.sorted.RG.bam \
 -O het/${SEX}/${FILENAME}.sorted.dedup.bam \
 -M het/${SEX}/${FILENAME}_marked_dedup_metrics.txt

rm het/${SEX}/${FILENAME}.sorted.RG.bam # Remove heavy sorted bam with read groups

