####################################################################################
#  name: 03_het_invar_join.R
#  desc: Joins het sites from script 01 and joins them to the list of all sites >10 cov identified in script 02
#  date: sept 2025
####################################################################################

library(tidyverse)
library(ggpubr)

# wd
setwd('/Users/benalston/Library/CloudStorage/GoogleDrive-btalston1@sheffield.ac.uk/My Drive/PhD/2024_sex_chr_id/Y_SNPDENSITY/het')

# data ----
# genotypes
M_genotypes <- read.csv('data/M_genotypes_ref_agnostic.csv')
F_genotypes <- read.csv('data/F_genotypes_ref_agnostic.csv')

# all sites (per 100kb window)
M_all_windows <- read_tsv('data/M_windows_filtered.tsv')
F_all_windows <- read_tsv('data/F_windows_filtered.tsv')
all_sites_joined <- rbind(
  mutate(M_all_windows, sex = 'male'),
  mutate(F_all_windows, sex = 'female')) %>% 
  rename('scaffold' = `#CHROM`)

# add window and sex info
M_genos_windows <- M_genotypes %>%
  group_by(scaffold, sample) %>% 
  mutate(window = (position %/% 100000) * 100000,
         sex = 'male')

F_genos_windows <- F_genotypes %>%
  group_by(scaffold, sample) %>% 
  mutate(window = (position %/% 100000) * 100000,
         sex = 'female')

# join genotype dfs and calculate number of het sites per sample per window
het_joined <- rbind(M_genos_windows,
                    F_genos_windows) %>% 
  group_by(scaffold, window, sample, sex) %>% 
  filter(het == 1) %>% 
  summarise(n_het = n())

# join with persite info
het_invar_sites_joined <- het_joined %>% left_join(all_sites_joined, by = c('scaffold', 'window', 'sex'))

# write csv ----
write_csv(het_invar_sites_joined, 'data/het_invar_sites.csv')
  
