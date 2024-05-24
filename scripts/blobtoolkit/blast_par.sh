#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --mem=6G
#SBATCH -c 6
#SBATCH --time=96:00:00
#SBATCH -e reports/blast/error_blast.txt
#SBATCH -o reports/blast/output_blast.txt
#SBATCH --array 1-500%50 # max jobs per user is 50

#########################################################################################
#	Script Name: busco.sh
#	Description: Runs blast on an assembly split into 500 smaller fasta files using fasta_splitter.sh, each containing numbers 1-500 in their filenames.
#	Author:      Ben Alston
#	Date:        April 2024
#########################################################################################

# Uses BLAST+ V2.15.0

wd=/mnt/parscratch/users/bip23bta/ref_genomes
species=whitei
assembly=2
blastdb=/mnt/parscratch/users/bip23bta/ref_genomes/blast_dbs/blast_nt_db/nt


cd $wd/$species/03-QC/blast/${assembly}_primary

apptainer exec ~/blast_latest.sif \
        blastn -db ${blastdb} \
       -query split/${assembly}_primary.part-${SLURM_ARRAY_TASK_ID}.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -out blast_out/${SLURM_ARRAY_TASK_ID}_blast.out

# If this script fails with an esoteric error about language localisation it is likely because you have forgotten to make the output directory specificed in the -out argument (/blast_out in this case)
# this script also requires you manually to condense the 500 individual blast output files e.g:
# cat blast_out/* > $assembly_name.blast.out
