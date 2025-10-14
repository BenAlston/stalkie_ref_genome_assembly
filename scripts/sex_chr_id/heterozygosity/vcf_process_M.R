####################################################################################
#  name: 01_vcf_process_F
#  desc: takes a filtered multi sample vcf file and calls het/hom sites agnostic to the ref or allele
#        subsequent scripts use this to calculate per window het. Seperate scripts for M and F
#  date: sept 2025
####################################################################################

# packages
library(vcfR)
library(tidyverse)

# wd
setwd('/Users/benalston/Library/CloudStorage/GoogleDrive-btalston1@sheffield.ac.uk/My Drive/PhD/2024_sex_chr_id/Y_SNPDENSITY/het/')

# read in vcf file
vcf <- read.vcfR('data/M_dal_7_filtered.vcf.gz')

# extract genotype (gt) information from vcf and clean sample names
sample_gts <- as.data.frame(vcf@gt)
sample_gts <- sample_gts %>% select(-FORMAT)
names(sample_gts) <- sub("-A*.*","", colnames(sample_gts)) # clean sample names

sample_gts_long <- sample_gts %>% 
  mutate(scaffold = vcf@fix[,"CHROM"], # add snp position and scaffold info
         position = vcf@fix[,"POS"]) %>% 
  pivot_longer(c(matches("Sample_")), # make dataframe long
               names_to = 'sample',
               values_to = 'genotype') %>% 
  mutate(geno1 = substr(genotype, 1,1), # extract the two genotypes from each sample
         geno2 = substr(genotype, 3,3),
         position = as.numeric(position)) %>% 
  select(-genotype)

# assigning hom or het irrespective of reference
sample_gts_ref_agnostic <- sample_gts_long %>% 
  mutate(het = 
           case_when(    
             geno1 == '.' & geno2 == '.' ~ NA, # missing data = NA
             geno1 ==  geno2 ~ 0, # two alleles are the same = hom
             geno1 != geno2 ~ 1 # two alleles are different = het
           )
  ) %>% 
  select(-geno1,-geno2)

# write csv
write_csv(sample_gts_ref_agnostic, 'data/M_genotypes_ref_agnostic.csv')

