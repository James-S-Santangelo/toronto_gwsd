# R script to smooth and average both genetic maps from Olsen et al. 2021
# Fits Shape Constrained Cubic P-spline to smoothed genetic maps

###############
###SETUP ####
###############

# Load requird packaes
library(tidyverse)
library(patchwork)
library(scam)

# Load sites files for interpolation
sites_allChroms <- snakemake@input[['sites']] %>%
  purrr::map_dfr(., read_delim, delim='\t', col_names=c('chrom', 'pos'))

# Load file with linkage map for all chromosomes
allMaps <- read_csv(snakemake@input[['linkage_map']])

# Creat list with chromosome names
chroms <- c('Chr01_Occ', 'Chr01_Pall', 'Chr02_Occ', 'Chr02_Pall', 'Chr03_Occ', 
            'Chr03_Pall', 'Chr04_Occ', 'Chr04_Pall', 'Chr05_Occ', 'Chr05_Pall', 
            'Chr06_Occ', 'Chr06_Pall', 'Chr07_Occ', 'Chr07_Pall', 'Chr08_Occ', 
            'Chr08_Pall')

# Load functions used throughout script
snakemake@source('./functions.R')

#########################
###FIT SCAM MODELS ####
#########################

# Fit SCAM model separately for each chromosome
scamFits_allChroms <- chroms %>%
  purrr::map_dfr(., get_fits, markers_df = allMaps, sites_df = sites_allChroms)

# Find minimum predicted cM
min_cM <- scamFits_allChroms %>% 
  group_by(chrom) %>% 
  summarise(min_cM = min(preds))

# Force minimum predicted cM to 0 for SCAM fits
scamFits_allChroms_off <- scamFits_allChroms %>% 
  group_by(chrom) %>% 
  mutate(preds = case_when(min(preds) < 0 ~ preds + abs(min(preds)),
                           min(preds) > 0 ~ preds - abs(min(preds))))

allMaps_offset <- allMaps %>% 
  left_join(., min_cM, by = 'chrom') %>%
  group_by(chrom) %>%
  mutate(cM = case_when(min_cM < 0 ~ cM + abs(min_cM),
                        min_cM > 0 ~ cM - abs(min_cM)))

# List of plots. Each plot is fit of SCAM model to chromosome markers
scamFits_plotList <- chroms %>%
  purrr::map(., scam_plot, markers_df = allMaps_offset, scam_fits_df = scamFits_allChroms_off)

# Combine list of plots into one
scamFits_allChroms_plot <- wrap_plots(scamFits_plotList, 4, 4)
ggsave(filename = snakemake@output[['scamFits_plot']], plot = scamFits_allChroms_plot,
       device = 'pdf', dpi = 600, width = 20, height = 20, units = 'in')

# Write interpolated genetic map for all chromosome to disk
scamFits_allChroms_off %>% 
    rename('cM' = 'preds') %>% 
    dplyr::select(pos, chrom, cM) %>%
    mutate(cM = as.numeric(cM)) %>%
    write_delim(., snakemake@output[['genMap_interp']], delim = '\t')




