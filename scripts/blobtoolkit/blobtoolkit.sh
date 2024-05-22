#!/bin/bash
#SBATCH -c 32 # CPU core
#SBATCH --mem=80G # memory required, up to 256G on standard nodes
#SBATCH -t 72:00:00 # time limit for job (format: days-hours:minutes:seconds)
#SBATCH -e reports/error_blobtoolkit.%J.txt
#SBATCH -o reports/output_blobtoolkit.%J.txt

#########################################################################################
#	Script Name: blobtoolkit.sh
#	Description: Runs blobtoolkit and filters an assembly to retain only contigs with metazoan blast hits or no blast hits
#                   Then runs BUSCO on the output
#	Author:      Ben Alston
#	Date:        April 2024
#########################################################################################

# needed modules: blobtoolkit2/3.1.6 samtools/1.17 minimap2/2.22-r1101
# blobtoolkit itself is called by an apptainer image (due to issues installing through conda)

# Load conda env containing minimap
module load Anaconda3
source activate blobtoolkit_run

# define apptainer images
blobtools='apptainer exec /users/bip23bta/blobtoolkit_latest.sif blobtools'
samtools='apptainer exec /users/bip23bta/samtools_latest.sif samtools'

# Variables
species=whitei
assembly=2
wd=/mnt/parscratch/users/bip23bta/ref_genomes/${species}
blobdir=${species}_${assembly}_primary

# Input data
INPUT_ASSEMBLY=$wd/02-hifiasm/${assembly}_hifiasm_output/*_primary.fa
HIFI_READS=$wd/data/long_read/200437_${assembly}-Cell?/*.fastq.gz
BUSCO=$wd/02-hifiasm/${assembly}_hifiasm_output/BUSCO_out_${assembly}/run_diptera_odb10/full_table.tsv
BLAST=$wd/03-QC/blast/${assembly}_primary/${assembly}_primary.blast.out

# set wd
cd $wd/03-QC/blobtools

# make blobdir
$blobtools create \
       --fasta $wd/02-hifiasm/${assembly}_hifiasm_output/*_primary.fa \
       ./$blobdir

### minimap2 to generate and add coverage

minimap2 -ax map-hifi -t 16 $INPUT_ASSEMBLY \
         $HIFI_READS \
        | samtools sort -@16 -O BAM -o data_files/$(basename $INPUT_ASSEMBLY).reads.bam

# samtools coverage does not work with older versions of samtools
$samtools coverage data_files/$(basename $INPUT_ASSEMBLY).reads.bam > data_files/$(basename $INPUT_ASSEMBLY).coverage.txt


# Add coverage to the blobdir
$blobtools add --text data_files/${assembly}_primary.fa.coverage.txt --text-header --text-cols '#rname=identifier,meandepth=my_reads_cov' \
    ./${species}_${assembly}_primary

$blobtools add --key plot.y=my_reads_cov ./$blobdir

# add blast hits to the blobdir

$blobtools add \
    --hits $BLAST \
    --hits-cols 1=qseqid,2=staxids,3=bitscore,5=sseqid,10=qstart,11=qend,14=evalue \
    --taxrule bestsum \
    --taxdump /mnt/parscratch/users/bip23bta/ref_genomes/taxdump \
    --replace \
    ./$blobdir

# generate fasta to be filtered
cat $INPUT_ASSEMBLY \
 > filtered_fastas/${assembly}_primary.fa

# filter the fasta
# keys=keep only the following. In this case, this keeps contigs with metazoan blast matches or no blast matches.
$blobtools filter \
    --param bestsum_kingdom--Keys=no-hit,Metazoa \
    --fasta filtered_fastas/${assembly}_primary.fa \
    --invert \
    $blobdir










