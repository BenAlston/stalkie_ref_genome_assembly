#!/bin/bash
#SBATCH --job-name=purge_dups
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH --time=24:00:00
#SBATCH -e reports/error_purge_dups_man_all.txt
#SBATCH -o reports/output_purge_dups_man_all.txt

#########################################################################################
#	Script Name: purge_dups_all.sh
#	Description: Runs purge dups on an assembly. Requires imput .fa, and the reads used to 
# 			       Generate the assembly
#	Author:      Ben Alston
#	Date:        May 2024
#########################################################################################

# purge_dups v1.2.5, minimap2 v2.27-r1193

# load packages
module load Anaconda3/2022.05
source activate purge_dups
bin=/users/bip23bta/purge_dups/bin # specifies location of pure dups scripts
minimap2='apptainer exec /users/bip23bta/minimap2_latest.sif minimap2' # using an apptainer image to call minimapup

# sample variables
SPECIES=whitei
ASSEMBLY=1

#filepaths
wd=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/03-QC/purge_dups/${ASSEMBLY}_primary
READS=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/data/long_read/200437_${ASSEMBLY}-Cell?/*.fastq.gz
REF_1=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/02-hifiasm/${ASSEMBLY}_hifiasm_output/${ASSEMBLY}_primary.fa

######## ---------- Round 1 ---------- ########
mkdir $wd
cd $wd
mkdir round_1
cd round_1

# map input reads to assembly
for i in $READS
do
  	$minimap2 -xasm20 $REF_1 $i | gzip -c - > $(basename $i).paf.gz
done

# generate cutoffs
${bin}/pbcstat *.paf.gz # (produces PB.base.cov and PB.stat files)
${bin}/calcuts PB.stat > cutoffs 2>calcults.log

# split assembly and self align:
${bin}/split_fa $REF_1 > $(basename $REF_1).split
minimap2 -xasm5 -DP $(basename $REF_1).split $(basename $REF_1).split | gzip -c - > $(basename $REF_1).split.self.paf.gz

# purge haplotigs
${bin}/purge_dups -2 -T cutoffs -c PB.base.cov $(basename $REF_1).split.self.paf.gz > dups.bed 2> purge_dups.log
${bin}/get_seqs -e dups.bed $REF_1

# merge hap and purged assemblies
cat hap.fa purged.fa > ${ASSEMBLY}_round-1_purged.fa


######## ---------- Round 2 ---------- ########

# new filepath
REF_2=/mnt/parscratch/users/bip23bta/ref_genomes/${SPECIES}/03-QC/purge_dups/${ASSEMBLY}_primary/round_1/${ASSEMBLY}_round-1_purged.fa

# change wd
cd $wd
mkdir round_2
cd round_2

# map input reads to round 1 purged assembly
for i in $READS
do
  	$minimap2 -xasm20 $REF_2 $i | gzip -c - > $(basename $i).paf.gz
done

# generate cutoffs 
${bin}/pbcstat *.paf.gz # (produces PB.base.cov and PB.stat files)
${bin}/calcuts PB.stat > cutoffs 2>calcults.log

# split assembly and self align:
${bin}/split_fa $REF_2 > $(basename $REF_2).split
minimap2 -xasm5 -DP $(basename $REF_2).split $(basename $REF_2).split | gzip -c - > $(basename $REF_2).split.self.paf.gz

# filter haplotigs
${bin}/purge_dups -2 -T cutoffs -c PB.base.cov $(basename $REF_2).split.self.paf.gz > dups.bed 2> purge_dups.log
${bin}/get_seqs -e dups.bed $REF_2




