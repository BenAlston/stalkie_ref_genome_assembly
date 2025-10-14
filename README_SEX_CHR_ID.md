# **Sex Chromosome Identification**
* Differences in coverage and heterozygosity can be used to identify the XY chromosomes
* Illumina short read data for the three species, 5 individuals per sex 
* this pipeline will be initially run on dalmanni
* wilkinson et al 2023 reported the tdal X to be 97.2 Mbp

## **Trimming & QC**
QC has already been done on the samples. Details of QC used by liverpool:
_The raw Fastq files are trimmed for the presence of Illumina adapter sequences using Cutadapt version 1.2.1. The option -O 3 was used, so the 3' end of any reads which match the adapter sequence for 3 bp. or more are trimmed._

_The reads are further trimmed using Sickle version 1.200 with a minimum window quality score of 20. Reads shorter than 15 bp. after trimming were removed. If only one of a read pair passed this filter, it is included in the R0 file. The output files from Cutadapt and Sickle are available here._

* Redoing qc so read length cuttoff can be 50
* Ran trimmomatic ([trim.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/QC/trim.sh)) with read cuttoff of 50
* [Multiqc report](https://github.com/BenAlston/stalkie_ref_genome_assembly/tree/main/lab_book/Data/multiqc_reports)
* Coverage and patterns of divergence suggest that samples A10 and A11 are potentially contaminated with eachother. These have been binned.

## **X chromosome Identification (MF Coverage Ratio)**
* Mapped to the female reference with bwa [01_mapping_bwa.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/coverage/01_mapping_bwa.sh). Stringent parameters to remove multi-mapping and retain only primary alignments.
* Generated a bed file with 100kb windows [02_bedtools_windows.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/coverage/02_bedtools_windows.sh)
* Use the resultant bam and bed files to extract coverage with for each sample with bedtools [03_bam_2_cov.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/coverage/03_bam_2_cov.sh)
* A sample ID col was then added to each of the resultant coverage TSV (tab separated) files manually in bash. These files were then collated and read into R:

~~~bash
# add a sample name column to the end of each coverage file:
for Sample in *_cov.tsv
do
echo $Sample
name=$(basename -s _cov.tsv $Sample)
awk -F'\t' -v OFS='\t' -v sample="$name" '{print $0, sample}' $Sample > ${name}_labeled.cov.tsv
done

# collate all files:
cat *_labeled.cov.tsv >> dal_7_cov_raw.tsv
~~~

### **X-autosomal cutoff**
* identify the X linked peak, then set an arbirary cutoff, 100kb windows within this cutoff will be considered X-linked. Rscript: [peak_id_final.R](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/coverage/peak_id_final.R)

stringency table:
| cutoff | num of contigs | total size (mb) | median size (mb) |
| ------ | -------------- | --------------- | ---------------- |
| \-0.1  | 45             | 2.93            | 0.03             |
| \-0.15 | 59             | 3.7             | 0.06             |
| \-0.2  | 75             | 4.34            | 0.03             |

renaming an assembly fasta:
~~~~bash
# renaming an assembly:
# the newnames.tsv file is a tab sep file with the format "old name"	"new name", (no header)

awk -F'\t' '{print "s/\\b" $1 "\\b/" $2 "/g"}' dal_7_newnames.tsv > rename.sed
sed -f rename.sed dal_7_scaffolded.fa > dal_7_renamed.fa
~~~~

# **Y identification: Degenerate Region**
* When mapping M and F reads to the M ref, regions where male read map only will be Y linked. Using the same scripts and parameters as the X
* Mapped with high strincency, same as X coverage id, same [coverage](https://github.com/BenAlston/stalkie_ref_genome_assembly/tree/main/scripts/sex_chr_id/coverage) scripts are used

# **Y identification: Diverged Region**
* Map M & F reads to M ref with less stringend parameters and call snps with [these scripts](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/heterozygosity)
