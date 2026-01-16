###########################################################################
# name: cov_process.R
# date: Nov 2025
# desc: Takes per contig coverage, calculates m:f cov ratio and plots
###########################################################################

# Packages ----
library(tidyverse)
library(readr)

# set wd
setwd('/Users/benalston/Library/CloudStorage/GoogleDrive-btalston1@sheffield.ac.uk/My Drive/PhD/2024_sex_chr_id/Y_ID/Y_cov_dal_6/cov_per_contig_samtools_new_params')

# ---- read in raw data ----
# sampleinfo
sampleinfo <- read_csv('data/sampleinfo.csv')

# read coverage and modify sample ids to match sampleinfo
cov_all_raw <- read_tsv('data/all_cov.tsv', 
                        col_names = c('contig','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq', 'id'))

# clean sample id and add sex
cov_all <- cov_all_raw %>% mutate(id = sub('-A.*', '', id)) %>%
  left_join(sampleinfo, by = c('id' = 'sample')) %>% rename(length = endpos)

#write.csv(cov_all, 'data/cov_processed.csv')

# blast info to sanity check
blast_raw <- read_csv('Y_dal_6_wgs_blast_pcr.csv')
blast_y_contigs <- blast_raw %>% filter(pident >= 99 & evalue == 0) %>%  select(sseqid, qseqid, pcr_validated) %>% 
  distinct()
blast_y_contigs_pcr <- blast_y_contigs %>% filter(pcr_validated == 'pcr_validated')

# ---- m:f cov ratio ----
# get mean cov per sex, convert na values to 0 coverage 
cov_ratio <- cov_all %>% group_by(sex, contig, length) %>% 
  summarise(mean_cov = mean(meandepth)) %>% 
  pivot_wider(names_from = sex, values_from=mean_cov) %>% 
  rename(mean_cov_m = Male,
         mean_cov_f = Female) %>% 
  mutate(mean_cov_m = ifelse(is.na(mean_cov_m), 0, mean_cov_m),
         mean_cov_f = ifelse(is.na(mean_cov_f), 0, mean_cov_f)) %>%
  mutate(mf_cov_ratio = (mean_cov_m+0.1)/(mean_cov_f+0.1))

cov_ratio_long <- cov_all %>% left_join(select(cov_ratio, contig,  mf_cov_ratio), by = 'contig')

# 15 of the 23 blast matched contigs are in my assembly, of these, 13 have been pcr validated.
cov_ratio %>% filter(contig %in% blast_y_contigs$sseqid) %>% nrow()
cov_ratio %>% filter(contig %in% blast_y_contigs_pcr$sseqid) %>% nrow()

# Y filters ----
ggplot(cov_ratio_long, aes(x = meandepth, fill = sex)) +
  geom_density(alpha = 0.5)
# no contigs with 0 f cov and >0 m cov

# get contigs where m and f contig overlaps
overlap_contigs <- cov_ratio_long %>% group_by(contig, sex) %>% 
  summarise(min = min(meandepth),
            max = max(meandepth)) %>% 
  pivot_wider(names_from = sex, values_from = c('min','max')) %>% 
  filter(max_Female >= min_Male) %>% pull(contig) %>% unique()

cov_ratio_long_filt <- cov_ratio_long %>% 
  filter(!contig %in% overlap_contigs) %>% 
  filter(mf_cov_ratio >= 2) %>%  # select contigs with cov ratio >= 2
  group_by(contig) %>% filter(any(meandepth >= 3)) # select contigs with highest cov >= 3

cov_ratio_long_filt %>% pull(contig) %>% unique() %>% length() # results in 20 contigs, 10mb
cov_ratio_long_filt %>% group_by(contig) %>% pull(length) %>% unique() %>% sum()

# only 2 of these have blast hits, both of which are pcr validated
cov_ratio_long_filt %>% filter(contig %in% blast_y_contigs$sseqid) %>% pull(contig) %>% unique() %>% length
cov_ratio_long_filt %>% filter(contig %in% blast_y_contigs_pcr$sseqid) %>% pull(contig) %>% unique() %>% length

cov_ratio_long_filt %>% group_by(sex) %>% summarise(meancov = mean(meandepth),
                                                    sd = sd(meandepth))



# long plots
cov_ratio_long_filt 

#pdf(file = 'cov_plot_long_filt.pdf', width = 14, height = 50)
cov_ratio_long_filt %>% 
  ggplot(aes(x = meandepth, y = fct_reorder(contig, mf_cov_ratio), fill = sex)) +
  geom_point(size = 0.5, alpha = 0.5) +
  stat_summary(fun.data = mean_se,
               geom = 'errorbar',
               width = 0.4) +
  stat_summary(geom = 'point',
               fun = mean,
               shape = 23,
               size = 3) +
  scale_fill_manual(values = c('red4', 'steelblue')) +
  theme_bw()
#dev.off()


cov_ratio_long %>% 
  filter(contig %in% unique(blast_y_contigs_pcr$sseqid)) %>% 
  ggplot(aes(x = meandepth, y = fct_reorder(contig, mf_cov_ratio), fill = sex)) +
  geom_point(size = 0.5, alpha = 0.5) +
  stat_summary(fun.data = mean_se,
               geom = 'errorbar',
               width = 0.4) +
  stat_summary(geom = 'point',
               fun = mean,
               shape = 23,
               size = 3) +
  scale_fill_manual(values = c('red4', 'steelblue')) +
  theme_bw()

