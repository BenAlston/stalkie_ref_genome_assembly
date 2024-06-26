# **Sex Chromosome Identification**
* Differences in coverage and heterozygosity can be used to identify the XY chromosomes
* Illumina short read data for the three species, 5 individuals per sex
* this pipeline will be initially run on whitei only

## **Trimming & QC**
QC has already been done on the samples. Details of QC used by liverpool:
_The raw Fastq files are trimmed for the presence of Illumina adapter sequences using Cutadapt version 1.2.1. The option -O 3 was used, so the 3' end of any reads which match the adapter sequence for 3 bp. or more are trimmed._

_The reads are further trimmed using Sickle version 1.200 with a minimum window quality score of 20. Reads shorter than 15 bp. after trimming were removed. If only one of a read pair passed this filter, it is included in the R0 file. The output files from Cutadapt and Sickle are available here._

* Redoing qc so read length cuttoff can be 50
* Ran trimmomatic ([trim.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/QC/trim.sh)) with read cuttoff of 50
* [Multiqc report](https://github.com/BenAlston/stalkie_ref_genome_assembly/tree/main/lab_book/Data/multiqc_reports)


~~~
female:A01,A02,A03,A04,A05
male:A06,A07,A08, A09, A10


~~~


**# Y chromosome Identification**



### **.1 mapping with bowtie2**
* Paired end reads: two .fastq.gz files per sample (R1 & R2) plus R0, an small file containing unpaired reads, disgarded.
* Run bowtie2 with --mp 10000 (effectivley removes missmatches by setting a high penalty for them)
* Ran [bowtie2_indexer.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/bowtie2_index.sh) on the female ref
* Ran [bowtie2_alignment](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/bowtie2_alignment.sh)
* mapped each of the individuals seperately to the female ref
 
### **2. extract per site coverage**
* [cov_calc.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/cov_calc.sh) (ongoing)
* Averaged across 5 kb windows using bedtools multicov
* generated merged .bed file with [bed_merge.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/bed_merge.sh)

### **3. Coverage ratio**
* In R, calculate log ratio of female to male coverage (per window) with dplyr
* used this to identify the X-linked reads

# **Y identification: Degenerate Region**
* WHen mapping M and F reads to the M ref, regions where male read map only will be Y linked. This is done in a similar way to the F reads

### **.1 mapping with bowtie2**
* Ran [bowtie2_indexer.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/bowtie2_index.sh) on male ref
* Ran [bowtie2_alignment](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/bowtie2_alignment.sh)
* mapped each of the individuals seperately to the male ref
 
### **2. extract per site coverage**
* [cov_calc.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/cov_calc.sh) 
* Calculated coverage across the genome in 5kb windows in exactly the same manner as the female ref alignment

### **3. Coverage ratio**
* In R, calculate log ratio of female to male coverage (per window)
* used this to identify the Y-linked reads

# **Y identification: Pseudoautosomal region**
