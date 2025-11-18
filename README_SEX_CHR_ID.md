# **Sex Chr ID**
intitially run on dalmanni only

## **Raw data preprocessing**
* 5M 5F per species
* Trimmed with trimmomatic [trim.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/QC/trim.sh)


## **X ID**
* ran [coverage scripts](https://github.com/BenAlston/stalkie_ref_genome_assembly/tree/main/scripts/sex_chr_id/coverage) to calculate mean coverage per contig and classify each unplaced contig in as X, autosomal or unknown
* Used to modify contig names

## **Y ID**
### **Coverage**
* ran [coverage scripts](https://github.com/BenAlston/stalkie_ref_genome_assembly/tree/main/scripts/sex_chr_id/coverage) on the male ref to rule out Y and autosomal linked contigs

### **Het**
* [SNP calling with GATK](https://github.com/BenAlston/stalkie_ref_genome_assembly/tree/main/scripts/sex_chr_id/heterozygosity/snp_calling)
1. SNPS are called with gatk, joint calling male and female samples seperatley. MAF = 0.1, nmissing = 0.9, mindepth = 10, maxdepth = 60 (depth is per sample per site)

2. From resultant vcfs, every het site is counted
3. coverage is calculated per site per sample and filtered (only sites with 10 cov across all samples are retained)
4. number of het sites and number of total sites are used to work out het per contig

* Now need to calculate the number of variant sites :
The problem:
- need a file containing number of het sites and number of called sites (inc invariant ones). In windows and per contig
- variant sites should be 'ref agnostic'
- currently extract all variant sites in R (sites where each allele is different). Then divide this by the total number of sites with >10 coverage in that window (same filter for snp calling)
- yaccine suggests just getting the raw number of het sites per window per individual, could factor in the total number of covered sites per ind if thats needed


  - keep the R script?
  - find a better way of calculating snp density. Maybe endit the vcf through R, then calc snpdensity with vcftools, if it can be determined that all called sites would be used in this calculation
