#!/bin/bash
#SBATCH --job-name=vcf_filter
#SBATCH --time=12:00:00
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -o reports/6_vcf_filtered_f.txt

#########################################################################################
#       Script Name: 06_vcf_filter.sh
#       Description: Filters vcf files with bcftools filter, see below for parameters
#       Author:      Ben Alston
#       Date:        Jun 2025
#########################################################################################

# ---- bcftools/1.20 and vcftools/0.1.16---- #
bcftools="apptainer exec /mnt/parscratch/users/bip23bta/docker_images/bcftools_1.20--h8b25389_1.sif bcftools"
vcftools="apptainer exec /mnt/parscratch/users/bip23bta/docker_images/vcftools_0.1.16--pl5321hdcf5f25_11.sif vcftools"

module load Anaconda3
source activate vcf_edit

# ---- Filepaths ---- #
VCF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/04_genotype_vcfs/${1}_dal_7_YID.vcf.gz
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/06_filter_vcf/${1}


mkdir -p $WD 
cd $WD 
#
# set filters
MAF=0.1
MISS=0.9
QUAL=20
MIN_DEPTH=10
MAX_DEPTH=60


#  filtering with vcftools
$vcftools \
        --gzvcf $VCF \
        --remove-indels \
        --maf $MAF \
        --minQ $QUAL \
        --minDP $MIN_DEPTH \
        --maxDP $MAX_DEPTH \
        --recode \
        --out ${1}_dal_7_filtered

$bcftools view -Oz -o ${1}_dal_7_filtered.vcf.gz ${1}_dal_7_filtered.recode.vcf
rm dal_7_filtered.vcf.gz
