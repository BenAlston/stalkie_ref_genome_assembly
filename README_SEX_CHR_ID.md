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
~~~bash
# append file id to coverage sheet
for file in *_cov.tsv; do
    id=$(basename "$file" _cov.tsv)
    awk -F'\t' -v OFS='\t' -v id="$id" '{print $0, id}' "$file" >> all_cov.tsv
done
~~~
* Ran [y_id.R](https://github.com/BenAlston/stalkie_ref_genome_assembly/blob/main/scripts/sex_chr_id/coverage/y_id.R) on the output
  - M:F Cov ratio > 2 & any cov > 3 gives 33mb across 114 contigs - can probably afford to make these more stringent
  - 15 of the 23 Y candidates from  Mahajan & Bachtrog (2017) are present in my male assembly. 12/15 are captured by my cov parameters. 11 of these are pcr validated (out of 14 total pcr validated Y sequences).
 

Refs: 
Mahajan, S. and Bachtrog, D., 2017. Convergent evolution of Y chromosome gene content in flies. Nature communications, 8(1), p.785.
