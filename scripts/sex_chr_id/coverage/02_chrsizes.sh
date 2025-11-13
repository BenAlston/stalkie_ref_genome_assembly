# file run with bash from command line
samtools='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/samtools_1.21--h50ea8bc_0.sif samtools'


# filepaths ----
REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/refs/dal_7_scaffolded.fa
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/05-sex_chr_id/01_X_ID/cov_per_contig/${1}

mkdir -p $WD 
cd $WD 

$samtools faidx $REF
cut -f1,2 ${REF}.fai | awk '{print $1":1-"$2}' > chrsizes.list
