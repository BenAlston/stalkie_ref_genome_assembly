#!/bin/bash
#SBATCH --job-name=vcf_qc
#SBATCH --time=24:00:00
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -o reports/5_vcf_qc.txt

#########################################################################################
#       Script Name: 05_vcf_qc.sh
#       Description: Filters vcf files with bcftools filter, see below for parameters
#       Author:      Ben Alston
#       Date:        Jun 2025
#########################################################################################

# this script is light enough that it can be run from command line, subsetting takes ~10mins, the rest is speedy

# ---- bcftools/1.20 and vcftools/0.1.16---- #
bcftools="apptainer exec /mnt/parscratch/users/bip23bta/docker_images/bcftools_1.20--h8b25389_1.sif bcftools"
vcftools="apptainer exec /mnt/parscratch/users/bip23bta/docker_images/vcftools_0.1.16--pl5321hdcf5f25_11.sif vcftools"

module load Anaconda3
source activate vcf_edit

# ---- Filepaths ---- #
VCF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/04_genotype_vcfs/${1}_dal_7_YID.vcf.gz
OUTPUT_VCF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/05_vcf_qc/${1}_dal_7_subset.vcf
SUBSET_VCF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/05_vcf_qc/${1}_dal_7_subset.vcf.gz
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/05_vcf_qc
OUT=${1}_dal_7_y_samples_subset


mkdir -p $WD 
cd $WD 

# subset vcf
$bcftools view $VCF | vcfrandomsample -r 0.012 > $OUTPUT_VCF
# compress vcf
bgzip $OUTPUT_VCF
# index vcf
$bcftools index $SUBSET_VCF

# ---- generate stats with vcftools ---- #
# Calculate allele frequency
$vcftools --gzvcf $SUBSET_VCF --freq2 --out $OUT --max-alleles 2
# Calculate mean depth per individual
$vcftools --gzvcf $SUBSET_VCF --depth --out $OUT
# Calculate mean depth per site
$vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT
# Calculate site quality
$vcftools --gzvcf $SUBSET_VCF --site-quality --out $OUT
# Calculate proportion of missing data per individual
$vcftools --gzvcf $SUBSET_VCF --missing-indv --out $OUT
# Calculate proportion of missing data per site
$vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT
# Calculate het
$vcftools --gzvcf $SUBSET_VCF --het --out $OUT
