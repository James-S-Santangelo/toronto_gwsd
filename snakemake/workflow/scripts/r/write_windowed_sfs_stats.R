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

write_delim(sfs_stats_windowed_df, snakemake@output[["sfs_df"]], delim = "\t")
