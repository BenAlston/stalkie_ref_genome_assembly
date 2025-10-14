library(tidyverse)
library(ggpubr)

setwd('/Users/benalston/Library/CloudStorage/GoogleDrive-btalston1@sheffield.ac.uk/My Drive/PhD/2024_sex_chr_id/Y_SNPDENSITY/het')
# read data
het_sites_raw <- read.csv('data/het_invar_sites.csv')

# calculate mean number of het sites per covered sites (depth > 10) in each window
het_sites_mean <- het_sites_raw %>% 
  mutate(het_proportion = n_het/sites) %>% 
  group_by(sex, window, scaffold) %>% 
  summarise(mean_het = mean(het_proportion)) 

# plot to check
het_sites_mean %>% 
  filter(scaffold == 'scaffold_2') %>% 
  ggplot(aes(x = window, y = mean_het, col = sex)) +
  geom_point()


# csum plots:
# unplaced contigsw ----
cumsum_vals <- het_sites_mean %>%
  group_by(scaffold) %>% 
  summarise(max_bp = max(window)) %>%
  arrange(desc(max_bp)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0))

het_csum <- het_sites_mean %>% 
  inner_join(cumsum_vals, by = 'scaffold') %>% 
  mutate(bp_cum = window + bp_add)

# sanity check plot
# just placed
het_csum_placed <- het_csum %>%  filter(scaffold %in% c(paste0('scaffold_',1:4)))
het_csum_x <-  het_csum %>%  filter(scaffold %in% c(paste0('scaffold_',2)))

ggplot(het_csum_placed, aes(x = bp_cum, y = mean_het, col = sex)) +
  geom_point(alpha = 0.8) +
  geom_vline(data = het_csum_placed %>% group_by(scaffold) %>% 
               summarise(scaff_end = max(bp_cum)),
             aes(xintercept = scaff_end)) +
  scale_color_manual(values = c('#eb4034','#02ab9a')) +
  theme(legend.position = 'none') +
  labs(title = 'proportion of het sites, 100kb windows', x = 'position (bp)') +
  theme_bw()

# all contigs and X plots
plot_grid()

ggarrange(
  ggplot(het_csum_placed, aes(x = bp_cum, y = mean_het, col = sex)) +
    geom_point(size = 0.5, alpha = 0.8) +
    geom_vline(data = het_csum_placed %>% group_by(scaffold) %>% 
                 summarise(scaff_end = max(bp_cum)),
               aes(xintercept = scaff_end)) +
    scale_color_manual(values = c('#eb4034','#02ab9a')) +
    theme(legend.position = 'none') +
    labs(title = 'proportion of het sites, all contigs 100kb windows', x = 'position (bp)') +
    theme_bw()
  ,
  ggplot(het_csum_x, aes(x = bp_cum, y = mean_het, col = sex)) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = c('#eb4034','#02ab9a')) +
    theme(legend.position = 'none') +
    labs(title = 'proportion of het sites, chr_X, 100kb windows', x = 'position (bp)') +
    theme_bw()
  ,
  nrow = 2,
  common.legend = TRUE
)





ggarrange(
ggplot(het_csum_placed, aes(mean_het, fill = sex)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c('#eb4034','#02ab9a')) +
  xlim(0, 0.03) +
#  ylim(0, 1800) +
  labs(title = 'mean het all placed scaffolds') +
  theme_bw()
,
ggplot(het_csum_x, aes(mean_het, fill = sex)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c('#eb4034','#02ab9a')) +
  xlim(0, 0.03) +
#  ylim(0, 1800) +
  labs(title = 'mean het chr_X') +
  theme_bw()
,
common.legend = TRUE
)








