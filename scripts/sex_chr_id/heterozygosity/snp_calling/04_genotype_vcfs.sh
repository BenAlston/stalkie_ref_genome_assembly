#!/bin/bash
#BATCH --job-name=genotype_vcf
#SBATCH --mem=100G
#SBATCH -c 4
#SBATCH --time=24:00:00
#SBATCH -o reports/4_vcf_genotype_output.txt

#########################################################################################
#       Script Name: 04_genotype_vcfs.sh
#       Description: Runs gatk GenotypeGVCFs to genotype combined GVCF files
#       Author:      Ben Alston
#       Date:        May 2024
#########################################################################################

# ---- Load GATK (versions:GATK/4.6.2.0,HTSJDK/3.0.1,Picard/2.27.5) ---- #
gatk='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/gatk_latest.sif gatk'
# working dir
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/04_genotype_vcfs
GENEDB=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/03_genomicsdb/${1}_db
REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/refs/dal_6_blobtools_scaffolded.fa
INTERVALS=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/02_Y_ID/snpdensity/het_dal_6/dal_6_intervals.list

# ---- Run script ---- #
mkdir -p $WD 
cd $WD 

# ---- GENOTYPE GVCFS SEPARATELY PER INTERVALS FOR EACH SPECIES ---- #
$gatk GenotypeGVCFs \
        --java-options "-Xmx15g" \
        -R $REF \
        -V gendb://${GENEDB} \
        -O ${1}_dal_6_YID.vcf.gz \
        --all-sites \
        --intervals ${INTERVALS}



# Extremely high memory useage for GenotypeGVCFs (200gb+) is an issue in older versions of gatk
# switching from gatk 4.3.0.0 to 4.6.2.0 solved the issue
