#!/bin/bash
#BATCH --job-name=hic_align
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --time=48:00:00
#SBATCH -o reports/bam_align_test.txt

# this is the script used in the arima genomics mapping pipeline

# filepaths
REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/7_omni-c_hifiasm_output_final/dal_7_raw.fa
READS1=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/data/omni_c/Tdalmanni-Female_R1_001.fastq.gz
READS2=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/data/omni_c/Tdalmanni-Female_R2_001.fastq.gz
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/04-scaffolding/1_mapping

# arima scripts (pulled directly from their git repo) & misc
STATS=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/04-scaffolding/scripts/mapping_pipeline/get_stats.pl
BAM_COMBINER=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/04-scaffolding/scripts/mapping_pipeline/two_read_bam_combiner.pl
FILTER=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/04-scaffolding/scripts/mapping_pipeline/filter_five_end.pl
CPU=20
MAPQ_FILTER=10

# load conda env containing samtools/1.3.1(wont work with newer versions) picard/3.4.0 and bwa/0.7.17
module load Anaconda3
source activate hi-c_mapping

# make and set wd
mkdir -p $WD 
cd $WD 

echo "### Step 0: Index reference"
bwa index -a bwtsw -p $(basename -s .fa $REF) $REF


echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t $CPU $REF $READS1 | samtools view -@ $CPU -Sb - > $(basename -s .fa $REF).1.bam
echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t $CPU $REF $READS2 | samtools view -@ $CPU -Sb - > $(basename -s .fa $REF).2.bam

echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h $(basename -s .fa $REF).1.bam | perl $FILTER | samtools view -Sb - > $(basename -s .fa $REF).1.filtered.bam

echo "### Step 2.B: Filter 5' end (2st)"
samtools view -h $(basename -s .fa $REF).2.bam | perl $FILTER | samtools view -Sb - > $(basename -s .fa $REF).2.filtered.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $BAM_COMBINER $(basename -s .fa $REF).1.filtered.bam $(basename -s .fa $REF).2.filtered.bam samtools  $MAPQ_FILTER | samtools view -bS -t ${REF}.fai - | samtools sort -@ $CPU -o $(basename -s .fa $REF).filtered.bam -

echo "### Step 3.B: Add read group"
picard AddOrReplaceReadGroups \
       INPUT=$(basename -s .fa $REF).filtered.bam \
       OUTPUT=$(basename -s .fa $REF).filtered.RG.bam \
       ID=$(basename -s .fa $REF) \
       LB=$(basename -s .fa $REF) \
       SM='stalkie_scaffolding' \
       PL=ILLUMINA \
       PU=none


echo "### Step 4: Mark duplicates"
picard MarkDuplicates \
       INPUT=$(basename -s .fa $REF).filtered.RG.bam \
       OUTPUT=$(basename -s .fa $REF).dedup.bam \
       METRICS_FILE=$(basename -s .fa $REF).txt \
       TMP_DIR=./ \
       ASSUME_SORTED=TRUE \
       VALIDATION_STRINGENCY=LENIENT \
       REMOVE_DUPLICATES=TRUE

samtools index $(basename -s .fa $REF).dedup.bam

perl $STATS $(basename -s .fa $REF).dedup.bam > $(basename -s .fa $REF).dedup.bam.stats
