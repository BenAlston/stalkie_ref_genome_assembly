#!/bin/bash
#BATCH --job-name=run_haphic
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=48:00:00
#SBATCH -o reports/run_haphic.txt

#########################################################################################
#       Script Name: run_haphic.sh
#       Description: runs haphic 
#       Author:      Ben Alston
#       Date:        April 2025
#########################################################################################

# load samtools and samblaster
module load Anaconda3
source activate haphic

# run haphic
HAPHIC=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/HapHiC/haphic
REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/dalmanni_6/6_omni-c_phasedhaps.fa
BAM=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/dalmanni_6/tdal_hic_filtered.bam
HAP1_GFA=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/6_omni-c_hifiasm_output/dalmanni_6_omni-c.asm.hic.hap1.p_ctg.gfa
HAP2_GFA=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/6_omni-c_hifiasm_output/dalmanni_6_omni-c.asm.hic.hap2.p_ctg.gfa


cd /mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/dalmanni_6/

$HAPHIC pipeline $REF $BAM 4 --gfa "$HAP1_GFA,$HAP2_GFA"
# 3 = number of chrs
