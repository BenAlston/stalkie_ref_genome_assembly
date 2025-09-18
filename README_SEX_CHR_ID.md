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
Get fasta seq names and lengths (for plotting):
~~~bash
cat file.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > contigs.tsv
~~~

### **Low coverage cutoff**
* [coverage histogram for whitei](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/lab_book/Data/sex_chr/whitei_cov_hist.jpg)
* symetrical outlier peaks on either side of the main peak are caused by regions with very low coverage, usually 1-4 per individual, symetry is determined by if its m or f coverage.
* Currently, samples with coverage values below <4 are considered 0, this has remomved the peaks. Still need to decide on a less arbirary threshold

### **X-autosomal cutoff**
* I have a multi-modal density distribution that looks like two normal distributions stuck together. 
* I essentially want to determine the probability that any given observation belongs to either peak, so I can set a justifiable cutoff
* doing this is kind of beyond me. I can fit models but have no idea if they are apropriate.

1. [Beatriz & Bachtrog et al (2015)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002078) identify the autosomal peak, divide that by two, and select contigs in a +-0.1 range around this.
2. The findZX pipeline doesn't say how they do this bit

# **Y identification: Degenerate Region**
* When mapping M and F reads to the M ref, regions where male read map only will be Y linked. Using the same scripts and parameters as the X
* Contigs with Ordinary male coverage and 0/near-0 female coverage can be considered X linked

# **Y identification: Diverged Region**
* Extract SNP density per KB from female alignments, this can be done with VCFtools
* Need to run gatk on the F genome to get a vcf
