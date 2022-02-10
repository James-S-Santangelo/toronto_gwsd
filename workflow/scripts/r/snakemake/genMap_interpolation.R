# R script to smooth and average both genetic maps from Olsen et al. 2021
# Fits Shape Constrained Cubic P-spline to smoothed genetic maps

###############
#### SETUP ####
###############

# Load requird packaes
library(tidyverse)
library(patchwork)
library(slider)
library(scam)

# Load functions used throughout script
snakemake@source('./functions.R')

# Assign input files to variables
DG_markers <- load_marker_blast(snakemake@input[['DG_markers']]) %>%
    mutate(pop = 'DG')
SG_markers <- load_marker_blast(snakemake@input[['SG_markers']]) %>%
    mutate(pop = 'SG')

# Load both genetic maps from Olsen et al. 
DG_genMap <- load_genMap(snakemake@input[['DG_genMap']]) %>%
    mutate(pop = 'DG')
SG_genMap <- load_genMap(snakemake@input[['SG_genMap']]) %>%
    mutate(pop = 'SG')

# Load files with marker names
DG_names <- read.csv(snakemake@input[['DG_names']], col.names = c('qseqid', 'name')) %>%
    mutate(pop = 'DG')
SG_names <- read.csv(snakemake@input[['SG_names']], col.names = c('qseqid', 'name')) %>%
    mutate(pop = 'SG')

# Load sites files for interpolation
sites_allChroms <- snakemake@input[['sites']] %>%
    purrr::map_dfr(., read_delim, delim='\t', col_names=c('chrom', 'pos'))

#################################################
#### FILTERING MARKERS AND REMOVING OUTLIERS ####
#################################################

# Combine markers from both mapping populations
all_markers <- bind_rows(SG_markers, DG_markers)

# Combine genetic maps from both populations
all_maps <- bind_rows(SG_genMap, DG_genMap)

# Combine marker names
all_names <- bind_rows(SG_names, DG_names)

# Find markers with only a single 100% match to genome
markers_onlyOneExact <- all_markers %>%
	group_by(prev_id, pop) %>%
	mutate(is_exact_match = ifelse(pident == 100, 1, 0)) %>%
	summarise(n = n(),
			n_exact = sum(is_exact_match)) %>%
	filter(n_exact == 1)

# Filter markers to only include those with a single match 
markers_blast_exact <- all_markers %>%
	filter(pident == 100) %>%
	filter(prev_id %in% markers_onlyOneExact$prev_id) %>%
	dplyr::select(qseqid, prev_chrom, prev_pos, prev_id, curr_chrom, curr_pos, curr_id, pop) %>%
	left_join(., all_names, by = c('qseqid', 'pop')) %>%
	left_join(., all_maps, by = c('name', 'pop'))
 
# Remove markers where Linkage group from crosses doesn't match chromosome from BLASTing to genome
# Identify outlier markers using marker-based, sliding window, plus distance to nearest markers
z_cut <- 1  # Z-score cutoff for outliers
num_markers <- 6  # Number of markers to left and right of focal marker in window
marker_dist_cut <- 1e7  # Maximum distance to nearest markers to be included
markers_blast_exact_filt <- markers_blast_exact %>%
	filter(LG == curr_chrom) %>%
	dplyr::select(curr_chrom, curr_pos, pop, cM) %>%
	rename("chrom" = "curr_chrom",
		 "pos" = "curr_pos") %>%
	group_by(chrom, pop) %>%
	arrange(pos, .by_group = TRUE) %>%

	# Points with no markers within `marker_dist_cut` will be excluded
	mutate(prev_pos = lag(pos),
		 next_pos = lead(pos),
		 prev_pos_dist = pos - prev_pos,
		 next_pos_dist = next_pos - pos,
		 lone_marker = case_when(is.na(prev_pos_dist) & next_pos_dist > marker_dist_cut ~ 1,
								 is.na(next_pos_dist) & prev_pos_dist > marker_dist_cut ~ 1,
								 prev_pos_dist > marker_dist_cut & next_pos_dist > marker_dist_cut ~ 1,
								 TRUE ~ 0)) %>%

	# Exclude markers above or below `z_cut` z-scores in windows considering `num_markers` to left and right
	mutate(win_mean = slide_dbl(cM, mean, .before = num_markers, .after = num_markers),
		 win_sd = slide_dbl(cM, sd, .before = num_markers, .after = num_markers),
		 marker_z = (cM - win_mean) / win_sd,
		 exclude = case_when(marker_z > z_cut | marker_z < -z_cut ~ 1,
							 lone_marker == 1 ~ 1,
							 TRUE ~ 0)) %>%

	# Remove DG map for chromosome 11 since missing markers for > 50% of chromosome
	filter(!(chrom == 'CM019111.1' & pop == 'DG')) %>%

	dplyr::select(chrom, pos, pop, cM, exclude)

# Remove outliers
markers_blast_exact_filt_noOut <- markers_blast_exact_filt %>%
	filter(exclude != 1)

# Plot markers, colored by mapping population, with outliers marked as X. All chromosomes
chroms <- markers_blast_exact_filt %>% pull(chrom) %>% unique()
markerPlot_byPop_list <- chroms %>%
	purrr::map(., plot_markers_byPop, df = markers_blast_exact_filt)

# Combine an save plots
markers_allChroms_byPop <- wrap_plots(markerPlot_byPop_list, 4, 4)
ggsave(filename = snakemake@output[['markers_byPop_plot']], plot = markers_allChroms_byPop, 
	   device = 'pdf', dpi = 600, width = 20, height = 20, units = 'in')

#########################################
#### SMOOTH MAPS AND FIT SCAM MODELS ####
#########################################

# Get mean cM across both maps in windows
winMeans_df <- chroms %>%
    purrr::map_dfr(., calculate_windowed_means, 
                   window_size = 5000000,
                   step = 2500000,
                   df = markers_blast_exact_filt_noOut)

# Fit SCAM model separately for each chromosome
scamFits_allChroms <- chroms %>%
    purrr::map_dfr(., get_fits, markers_df = winMeans_df, sites_df = sites_allChroms) 

# Find minimum predicted cM
min_cM <- scamFits_allChroms %>% 
	group_by(chrom) %>% 
	summarise(min_cM = min(preds))

# Force minimum predicted cM to 0 for SCAM fits
scamFits_allChroms <- scamFits_allChroms %>% 
	group_by(chrom) %>% 
	mutate(preds = case_when(min(preds) < 0 ~ preds + abs(min(preds)),
							 min(preds) > 0 ~ preds - abs(min(preds))))

# Offset markers by same amount as fits. Only for plotting
winMeans_df <- winMeans_df %>% 
	left_join(., min_cM, by = 'chrom') %>% 
	group_by(chrom) %>% 
	mutate(win_mean_cM_off = case_when(min_cM < 0 ~ win_mean_cM + abs(min_cM),
									   min_cM > 0 ~ win_mean_cM - abs(min_cM)))

# List of plots. Each plot is fit of SCAM model to chromosome markers
scamFits_plotList <- chroms %>%
    purrr::map(., scam_plot, markers_df = winMeans_df, scam_fits_df = scamFits_allChroms)

# Combine list of plots into one
scamFits_allChroms_plot <- wrap_plots(scamFits_plotList, 4, 4)
ggsave(filename = snakemake@output[['scamFits_plot']], plot = scamFits_allChroms_plot,
       device = 'pdf', dpi = 600, width = 20, height = 20, units = 'in')

# Write interpolated genetic map for all chromosome to disk
scamFits_allChroms %>% 
    rename('cM' = 'preds') %>% 
    dplyr::select(pos, chrom, cM) %>%
    mutate(cM = as.numeric(cM)) %>%
    write_delim(., snakemake@output[['genMap_interp']], delim = '\t')




