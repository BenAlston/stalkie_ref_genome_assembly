#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --mem=6G
#SBATCH -c 6
#SBATCH --time=96:00:00
#SBATCH -e reports/blast/error_blast.txt
#SBATCH -o reports/blast/output_blast.txt
#SBATCH --array 1-500%50 # max jobs per user is 50

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
