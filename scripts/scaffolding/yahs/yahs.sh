#!/bin/bash
#SBATCH --job-name=yahs
#SBATCH --mem=100G
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH -o reports/yahs.txt



yahs='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/yahs_1.2.2--h577a1d6_1.sif yahs'
juicer='apptainer exec /mnt/parscratch/users/bip23bta/docker_images/yahs_1.2.2--h577a1d6_1.sif juicer'


REF=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/7_omni-c_hifiasm_output_final/dal_7_raw.fa
BAM=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/dal_7/tdal_hic_aligned.bam
WD=/mnt/parscratch/users/bip23bta/ref_genomes/dalmanni/02-hifiasm/haphic/dal_7/yahs

cd $WD 

#$yahs $REF $BAM

#($juicer pre *.bin yahs.out_scaffolds_final.agp ${REF}.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

cut -f1-2 ${REF}.fai > yahs.out_scaffolds_final.chrom.sizes

module load Anaconda3
source activate juicertools
#(juicer_tools pre alignments_sorted.txt out.hic.part yahs.out_scaffolds_final.chrom.sizes) && (mv out.hic.part out.hic)

#$juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp ${REF}.fai >out_JBAT.log 2>&1

#$juicer pre -a -o out_JBAT *.bin yahs.out_scaffolds_final.agp ${REF}.fai  > out_JBAT.log 2>&1
export _JAVA_OPTIONS="-Xmx100G"
(juicer_tools pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
