# installs and decompresses the blast nucleotide database, which takes up ~500gb
# run on command line, though would be a bit quicker as a slurm script
cd /mnt/parscratch/users/bip23bta/ref_genomes/blast_nt_db

# used an apptainer image because conda hates me
apptainer shell ~/blast_latest.sif 
update_blastdb.pl --decompress nt
