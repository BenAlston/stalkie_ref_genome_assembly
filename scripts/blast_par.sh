
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
assembly=1

cd $wd/$species/03-QC/blast/${assembly}_primary

apptainer exec ~/blast_latest.sif \
        blastn -db ${wd}/blast_nt_db/nt \
       -query split/${assembly}_primary.part-${SLURM_ARRAY_TASK_ID}.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -out blast_out/${SLURM_ARRAY_TASK_ID}_blast.out
