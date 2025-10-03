#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=24:00:00
#SBATCH -o reports/Y_run/%J_mapping.txt
#SBATCH --array 0-8%9

#########################################################################################
#       Script Name: 1_mapping.sh
#       Description: Aligns short reads to a ref, runs picard addorreplacereadgroups and markduplicates
#                    outputs a_sorted.dedup.bam file. Cov from resultant files will be used to id sex-linked contigs
#       Author:      Ben Alston
#       Date:        Sep 2025
#########################################################################################

# load bowtie2/2.4.1 and samtools/1.3.1
bwa='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/bwa_v0.7.17_cv1.sif bwa'
samtools='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/samtools_1.21--h50ea8bc_0.sif samtools'

# ---- Read array from file of filenames ---- #
FILEPREFIXES=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/sample_list.txt
readarray -t array < $FILEPREFIXES # a file containing a list of file prefixes 
FILENAME="${array[$SLURM_ARRAY_TASK_ID]}"

# ----filepaths ---- #
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/bams
READS_1=/mnt/parscratch/users/bip23bta/ref_genomes/*/data/short_read/Trimmed/${FILENAME}/*_L001_R1_001.fastq.gz
READS_2=/mnt/parscratch/users/bip23bta/ref_genomes/*/data/short_read/Trimmed/${FILENAME}/*_L001_R2_001.fastq.gz
REF_INDEX=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/refs/dal_6_blobtools_scaffolded.fa


#################
# bwa mem align #
#################
mkdir -p $WD 
cd $WD 

# reads are mapped strignently to eliminate multi mapping.
$bwa mem -t 4 -B 40 -O 60 -E 10 -L 100 $REF_INDEX $READS_1 $READS_2 | $samtools view -b -f 2 -F 2316 | $samtools sort -o ${FILENAME}_sorted.bam
# -- BWA parameters
# -B = missmatch penalty
# -O = Gap open penalty
# -E = Gap extension penalty
# -L = Clipping penalty
# Samtools parameters
# -F 2316 = Exclude reads with flags 2316 (Read unmapped, mate unmapped, not primary alignment, supplementary alignment)

###############################
# GATK addorreplacereadgroups #
###############################

module load GATK

RGID=$(zcat ${READS_1} | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4) # Extract RGID from read header
RGLB="Solexa"
RGSM="$FILENAME"
RGPU=$RGID.$RGLB
RGPL="Illumina"

gatk AddOrReplaceReadGroups \
        --I ${FILENAME}_sorted.bam \
        --O ${FILENAME}_sorted.RG.bam \
        --RGID $RGID \
        --RGLB $RGLB \
        --RGPL $RGPL \
        --RGPU $RGPU \
        --RGSM $RGSM \
        --SORT_ORDER coordinate \
        --CREATE_INDEX TRUE


#########################
# Picard MarkDuplicates #
#########################

mkdir -p TMP

gatk MarkDuplicates \
 --REMOVE_DUPLICATES TRUE \
 --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
 --VALIDATION_STRINGENCY SILENT \
 --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
 --TMP_DIR TMP \
 -I ${FILENAME}_sorted.RG.bam \
 -O ${FILENAME}_sorted.dedup.bam \
 -M ${FILENAME}_marked_dedup_metrics.txt

 ---- Index bam ---- #

$samtools index ${FILENAME}_sorted.dedup.bam

# ---- Remove unwanted Bams to save space ---- #
rm ${FILENAME}_sorted.RG.bam
rm ${FILENAME}_sorted.bam

