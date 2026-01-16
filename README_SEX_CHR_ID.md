# **Sex Chr ID**
intitially run on dalmanni only

## **Raw data preprocessing**
* 5M 5F per species
* Trimmed with trimmomatic [trim.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/QC/trim.sh)

## **Mapping**
* Mapping with bwa: [01_mapping_bwa.sh](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/coverage/01_mapping_bwa.sh).
  - Samflags 260, 0x2 (proper pair), -q 20
* This is run on both the female ref (X id) and male ref (Y id)

## **X ID**
* ran [coverage scripts](https://github.com/BenAlston/stalkie_ref_genome_assembly/tree/main/scripts/sex_chr_id/coverage) 02-04 to calculate mean coverage per contig and classify each unplaced contig in as X, autosomal or unknown
* Used to modify contig names

## **Y ID**
### **Coverage**
* ran [coverage scripts](https://github.com/BenAlston/stalkie_ref_genome_assembly/tree/main/scripts/sex_chr_id/coverage)
* Ran [y_id.R](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/coverage/y_id.R) on the output
