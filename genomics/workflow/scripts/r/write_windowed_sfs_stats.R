# Script to write windowed fst/theta and xpnsl dataframes

# Load required packages
library(tidyverse)

###################
#### SFS STATS ####
###################

# Function to load windowed Fst
load_windowed_fst <- function(path){
    colnames <- c("region", "Chr", "WinCenter", "nSites_fst", "fst")
    df <- suppressMessages(read_delim(path, delim = '\t', skip = 1, col_names = colnames)) %>%
    return(df)
}

# Load and clean windowed fst data
fst_df <- snakemake@input["fst"] %>%
    purrr::map_dfr(., load_windowed_fst)
names(fst_df) <- c("region", "Chr", "WinCenter", "nSites_fst", "fst")
fst_df_mod <- fst_df %>%
    mutate(fst = ifelse(fst < 0, 0, fst)) %>%
    mutate(region = str_extract(region, pattern = "(?<=\\()\\d+,\\d+(?=\\)$)")) %>%
    separate(region, into = c("start", "end"), sep = ",") %>%
    mutate(chrom_pos = paste0(Chr, "_", WinCenter))

# Load and clean urban and rural thetas
thetaU_df <- purrr::map_dfr(snakemake@input["thetaU"], read_delim, delim = "\t") %>%
    mutate(habitat = "Urban")
thetaR_df <- purrr::map_dfr(snakemake@input["thetaR"], read_delim, delim = "\t") %>%
    mutate(habitat = "Rural")
theta_df <- bind_rows(thetaU_df, thetaR_df) %>%
    rename(region = 1) %>%
    mutate(tp_scaled = tP / nSites,
           habitat = habitat) %>%
    dplyr::select(-tF, -tH, -tL, -fuf, -fud, -fayh, -zeng, -tW, -tP, -starts_with("#")) %>%
    mutate(region = str_extract(region, pattern = "(?<=\\()\\d+,\\d+(?=\\)$)")) %>%
    separate(region, into = c("start", "end"), sep = ",") %>%
    mutate(chrom_pos = paste0(Chr, "_", WinCenter))
thetas_df_wide <- theta_df %>%
    pivot_wider(values_from = c("Tajima", "nSites", "tp_scaled"), names_from = "habitat")

# Combine fst and theta dataframes
sfs_stats_windowed_df <- thetas_df_wide %>%
    filter(chrom_pos %in% fst_df_mod$chrom_pos) %>%
    left_join(., fst_df_mod, by = c("chrom_pos", "Chr", "start", "end", "WinCenter")) %>%
    mutate(delta_tp_ur = tp_scaled_Urban - tp_scaled_Rural,
           delta_td_ur = Tajima_Urban - Tajima_Rural)

# Identify outliers across genome
nSites_thresh <- as.numeric(snakemake@params[['nSites']])

win_sfs_df_filt <- sfs_stats_windowed_df %>%
    filter_at(vars(starts_with('nSites')), ~ . >= nSites_thresh)

sprintf("%s of %s Fst windows remaining", nrow(win_sfs_df_filt), nrow(sfs_stats_windowed_df))

fst_quant_filt <- quantile(win_sfs_df_filt %>% pull(fst), probs = c(0.99))
tp_quant_filt <- quantile(win_sfs_df_filt %>% pull(delta_tp_ur), probs = c(0.01, 0.99))
td_quant_filt <- quantile(win_sfs_df_filt %>% pull(delta_td_ur), probs = c(0.01, 0.99))

win_sfs_df_filt <- win_sfs_df_filt %>%
    mutate(fst_outlier = ifelse(fst >= fst_quant_filt, 1, 0),
           tp_outlier = ifelse(delta_tp_ur <= tp_quant_filt[1] | delta_tp_ur >= tp_quant_filt[2], 1, 0),
           td_outlier = ifelse(delta_td_ur <= td_quant_filt[1] | delta_td_ur >= td_quant_filt[2], 1, 0),
           all_outlier = ifelse(fst_outlier == 1 & tp_outlier == 1 & td_outlier == 1, 1, 0)) %>%
    dplyr::select(chrom_pos, Chr, start, end, WinCenter, fst, delta_tp_ur, delta_td_ur, contains('_outlier'))

# Add habitat under selection based on difference in pi and Tajima's D
win_sfs_df_filt <- win_sfs_df_filt %>% 
    mutate(direction = case_when(delta_tp_ur < 0 & delta_td_ur < 0 ~ 'Urban sel',
                                 delta_tp_ur > 0 & delta_td_ur > 0 ~ 'Rural sel',
                                 TRUE ~ 'Weird'))

write_delim(win_sfs_df_filt, snakemake@output[["sfs_df"]], delim = "\t")
