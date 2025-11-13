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
* Not confident in the above workflow, needs redoing
* Now need to calculate the number of variant sites 
