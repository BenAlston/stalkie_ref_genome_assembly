################################################################################################
#       Cov peak calc
#       Ben Alston
#       Sep 2025
################################################################################################

# packages
library(tidyverse)
library(readr)
library(cowplot)

# library
setwd( "/Users/benalston/Library/CloudStorage/GoogleDrive-btalston1@sheffield.ac.uk/My Drive/PhD/2024_sex_chr_id/X_cov_dal_7")

# read data
placed_contigs <- read.csv('data/mf_cov_ratio_placed_100kb.csv')
unplaced_contigs <- read.csv('data/mf_cov_ratio_unplaced_100kb.csv')


# plot density
plot_grid(
ggplot(placed_contigs, aes(x = mf_cov_ratio)) +
  geom_density() +
  labs(title = 'cov_ratio_placed') +
  theme_bw()
,
ggplot(unplaced_contigs, aes(x = mf_cov_ratio)) +
  geom_density() +
  labs(title = 'cov_ratio_unplaced') +
  theme_bw()
,
nrow = 1
)
# ------------------------- identify cov peaks ------------------------- # ----
# find max peak in density distribution (autosomal peak)

# --- input data 
input_data <- unplaced_contigs

# ---- first extract density dist
cov_density <- density(input_data$mf_cov_ratio)
xvals = cov_density$x
yvals = cov_density$y

# ---- get max value in y (density), find it's corresponding X cord: 
#
# diff(): calculates differences between consecutive elements in a numeric vec
# sign(): coverts these to 1 -1 or 0, depending on if they are +ve -ve or 0
# which(): identifies positions in a vector are true for a specific condition, in this case sign == -2
# sign will only = -2 if a positive value is next to a negative value
# in other words, the slope has changed from positive to negative
# we add 1 at the end because diff reduces the size of the vec
cov_peaks_index <- which(diff(sign(diff(cov_density$y))) == -2) + 1

# extract the cov ratio and density values that correspond to these peaks
peaks <- data.frame(cov_ratio = xvals[cov_peaks_index], density = yvals[cov_peaks_index])
# select the top 2 (which will always be the X and autosomal peaks in this case)
peaks <- peaks %>% arrange(desc(density)) %>% top_n(2) %>% pull(cov_ratio)
peaks
# plot to check
# arbirary 0.1 cutoff
cutoff = 0.1

ggplot(input_data, aes(x = mf_cov_ratio)) +
  geom_rect(aes(xmin = peaks[2]-cutoff, xmax = peaks[2]+cutoff, ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.1) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = peaks, col = 'red', linetype ='dashed') +
  xlim(0,2) +
  theme_bw()


# ------------------------- assigning X-linked contigs ------------------------- # ----

# define contigs within cutoffs as x or autosomal linked
unplaced_contigs_linkage <- unplaced_contigs %>% 
  mutate(linkage = case_when(
    between(mf_cov_ratio, peaks[2]-cutoff,peaks[2]+cutoff) ~ 'X',
    .default = "aut"
    ))

# calculate as percentage
unplaced_percentage <- unplaced_contigs_linkage %>% 
  group_by(contig) %>% 
  mutate(contig_length = length(contig)) %>% 
  group_by(contig, linkage) %>% 
  summarise(percent = length(linkage)/unique(contig_length)*100,
            length = unique(length)) %>% 
  arrange(linkage)

# plots ----
# plot larger contigs
ggplot(filter(unplaced_percentage, length >= 0.5e6), aes(x = fct_reorder(contig, linkage), y = percent, fill = linkage)) +
  geom_col() +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# plot all contigs
ggplot(unplaced_percentage, aes(x = fct_reorder(contig, linkage), y = percent, fill = linkage)) +
  geom_col() +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# plot coverage for single scaffold
scaffold = 'scaffold_10'
ggplot(filter(unplaced_contigs_linkage, contig == scaffold), aes(x = mf_cov_ratio, fill = contig)) +
  geom_density() +
  xlim(0,2) +
  labs(title = paste0(scaffold)) +
  theme_bw() +
  theme(legend.position = 'none')

# x_contig info ----

# extract 100% x-linked contigs
x_linked_unplaced_contigs <- unplaced_percentage %>% filter(linkage == 'X' & percent == 100)
x_linked_unplaced_contigs %>% pull(length) %>% sum()
x_linked_unplaced_contigs %>% summary()
#2,931,499bp in 45  potentially x linked unplaced contigs, median length = 29,690, max = 834,069
# worth considering this is < 1% of the genome
# get list of x-linked contigs
x_linked_contig_names <- x_linked_unplaced_contigs$contig






