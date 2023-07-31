# Script to write windowed fst/theta and xpnsl dataframes

# Load required packages
library(tidyverse)

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

#" Calculate mean cM of markers in window
#"
#" @param chrom_name Character vector with name of chromosome
#" @param window_size Size of window in bp
#" @param step Number of base pairs to shift window
#" @param df Dataframe with windowed markers
#"
#" @return Dataframe with mean cM in windows
calculate_windowed_stats <- function(df, window_size, step, thresh){
    chrom <- df %>% pull(Chr) %>% unique()
    winStarts <- seq(from = min(df$pos), to = max(df$pos), by = step)
    mat <- matrix(0, nrow = length(winStarts), ncol = 13)
    for(i in 1:length(winStarts)){
        start <- winStarts[i]
        end <- start + step
        df_filt <- df %>% filter(pos >= start & pos < end)
        winID <- i
        winCenter <- start + (step / 2)
        mean <- suppressWarnings(mean(df_filt$normxpnsl))
        max <- suppressWarnings(max(df_filt$normxpnsl))
        min <- suppressWarnings(min(df_filt$normxpnsl))
        n <- nrow(df_filt)
        num_gt_thresh <- sum(df_filt$normxpnsl > thresh)
        num_lt_thresh <- sum(df_filt$normxpnsl < -thresh)
        gt_frac <- sum(num_gt_thresh) / n
        lt_frac <- sum(num_lt_thresh) / n
        stats <- c(chrom, winID, start, end, winCenter, mean, max, min, n, num_gt_thresh, num_lt_thresh, gt_frac, lt_frac)
        mat[i, ] <- stats
    }
    stats_df <- as.data.frame(mat)
    names(stats_df) <- c("Chr", "winID", "start", "end", "winCenter", "mean", "max", "min", "n", "num_gt_thresh", "num_lt_thresh", "gt_frac", "lt_frac")
    return(stats_df)
}

load_xpnsl_norm <- function(path){
    chrom_name <- str_extract(basename(path), pattern = '.+(?=_Urban)')
    df <- suppressMessages(read_delim(path, delim = "\t")) %>%
        mutate(Chr = chrom_name)
    return(df)
}

# Load raw xpnsl data
xpnsl_norm_df <- snakemake@input[["xpnsl"]] %>%
    purrr::map_dfr(., load_xpnsl_norm) %>%
    rename("normxpnsl" = "normxpehh") %>%
    dplyr::select(-id)

# Calculate windowed xpnsl
window_size <- as.numeric(snakemake@params["winsize"])
step <- as.numeric(snakemake@params["winsize"])
thresh <- 2
xpnsl_windows <- xpnsl_norm_df %>%
    group_split(Chr) %>%
    purrr::map_dfr(., calculate_windowed_stats, window_size = window_size, step = step, thresh = thresh)
write_delim(xpnsl_windows, snakemake@output[["xpnsl_df"]], delim = "\t")
