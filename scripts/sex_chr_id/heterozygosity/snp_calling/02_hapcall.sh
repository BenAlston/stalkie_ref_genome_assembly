#!/bin/bash
#BATCH --job-name=varcall
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH --time=96:00:00
#SBATCH -o reports/Y_run/%J_02_hapcall.txt
#SBATCH --array 0-8%9

#########################################################################################
#       Script Name: 02_hapcall.sh
#       Description: variant calling with GATK requires sorted bam files, outputs vcf files
#       Author:      Ben Alston
#       Date:        Aug 2024
#########################################################################################

# ---- Read array from file of filenames ---- #
FILEPREFIXES=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/sample_list.txt
readarray -t array < $FILEPREFIXES # a file containing a list of file prefixes 
FILENAME="${array[$SLURM_ARRAY_TASK_ID]}"

# load load GATK/4.3.0
module load GATK
samtools='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/samtools_1.21--h50ea8bc_0.sif samtools'

# working dir
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/02_varcall

# reads and index
BAM=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/01_bams/${FILENAME}_sorted.dedup.bam
REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/refs/dal_6_blobtools_scaffolded.fa
OUTPUT=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/02_varcall

# ---- index ref ---- #

cd /mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/refs/

gatk CreateSequenceDictionary \
        --R $REF

$samtools faidx $REF
# ---- Run gatk haplotypecaller ---- # 
cd $WD 

mkdir -p $OUTPUT

gatk HaplotypeCaller \
-R $REF \
-I $BAM \
-O $OUTPUT/$(basename -s _sorted.dedup.bam ${BAM}).g.vcf.gz \
-ERC GVCF \
--dont-use-soft-clipped-bases \
--minimum-mapping-quality 20 \
--base-quality-score-threshold 20 \
--output-mode EMIT_ALL_CONFIDENT_SITES \
--native-pair-hmm-threads 12

# ---- reblock vcf to reduce memory usage downstream ---- #

gatk ReblockGVCF \
   -GQB 20 -GQB 30 -GQB 40 --floor-blocks \
   -R $REF \
   -V $OUTPUT/$(basename -s _sorted.dedup.bam ${BAM}).g.vcf.gz \
   -O $OUTPUT/$(basename -s _sorted.dedup.bam ${BAM}).g.rb.vcf.gz
