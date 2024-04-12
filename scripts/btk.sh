#!/bin/bash
#SBATCH -c 32 # CPU core
#SBATCH --mem=80G # memory required, up to 256G on standard nodes
#SBATCH -t 72:00:00 # time limit for job (format: days-hours:minutes:seconds)
#SBATCH -e reports/error_blobtoolkit.%J.txt
#SBATCH -o reports/output_blobtoolkit.%J.txt

# needed modules: blobtoolkit2/3.1.6 samtools/1.17 minimap2 blast
# blobtoolkit was installed with: <pip install "blobtoolkit[full]"> so im not sure its part of the conda env

module load Anaconda3
source activate blobtoolkit_run
blobtools='apptainer exec /users/bip23bta/blobtoolkit_latest.sif blobtools'
samtools='apptainer exec /users/bip23bta/samtools_latest.sif samtools'

# VARIABLES
species=whitei
assembly=2
wd=/mnt/parscratch/users/bip23bta/ref_genomes/${species}
INPUT_ASSEMBLY=${wd}/02-hifiasm/${assembly}_hifiasm_output/*_primary.fa
HIFI_READS=${wd}/data/long_read/200437_${assembly}-Cell?/*.fastq.gz
blobdir=${species}_${assembly}_primary
BUSCO=$wd/02-hifiasm/${assembly}_hifiasm_output/BUSCO_out_${assembly}/run_diptera_odb10/full_table.tsv

# set wd
cd ${wd}/03-QC/blobtools

# make blobdir
#apptainer exec ~/blobtoolkit_latest.sif blobtools create \
#	--fasta $INPUT_ASSEMBLY \
#	$blobdir

### minimap2 to generate and add coverage

#minimap2 -ax map-hifi -t 16 $INPUT_ASSEMBLY \
#         $HIFI_READS \
#        | samtools sort -@16 -O BAM -o data_files/$(basename $INPUT_ASSEMBLY).reads.bam

# samtools coverage does not work with older versions of samtools
#$samtools coverage data_files/$(basename $INPUT_ASSEMBLY).reads.bam > data_files/$(basename $INPUT_ASSEMBLY).coverage.txt

#$blobtools add --text data_files/$(basename $INPUT_ASSEMBLY).coverage.txt --text-header --text-cols '#rname=identifier,meandepth=my_reads_cov' ./$blobdir

#$blobtools add --key plot.y=my_reads_cov ./$blobdir

blobtools add \
    --hits /mnt/parscratch/users/bip23bta/ref_genomes/whitei/03-QC/blast/1_primary/1_blast.out \
    --taxrule bestsum \
    --taxdump /mnt/parscratch/users/bip23bta/ref_genomes/taxdump \
    --replace \
    whitei_1_primary

# generate fa to be filtered, not sure why they did this rather than cp
cat /mnt/parscratch/users/bip23bta/ref_genomes/whitei/02-hifiasm/1_hifiasm_output/1_primary.fa \
 > 1_primary.fa-uncont.fa

# keys=no-hit means exclude

blobtools filter \
    --param bestsum_superkingdom--Keys=Eukaryota \
    --invert \
    --fasta 1_primary.fa-uncont.fa \
    whitei_1_primary
