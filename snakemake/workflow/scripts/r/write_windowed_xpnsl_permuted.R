# Script to write windowed fst/theta and xpnsl dataframes

# Load required packages
library(tidyverse)

################
#### XP-nSL ####
################

calculate_windowed_stats <- function(df, window_size, step, thresh){
    chrom <- df %>% pull(Chr) %>% unique()
    perm <- df %>% pull(iter) %>% unique()
    winStarts <- seq(from = 0, to = max(df$pos) + window_size, by = step)
    mat <- matrix(0, nrow = length(winStarts), ncol = 14)
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
        stats <- c(chrom, perm, winID, start, end, winCenter, mean, max, min, n, num_gt_thresh, num_lt_thresh, gt_frac, lt_frac)
        mat[i, ] <- stats
    }
    stats_df <- as.data.frame(mat)
    names(stats_df) <- c("Chr", "iter", "winID", "start", "end", "winCenter", "mean", "max", "min", "n", "num_gt_thresh", "num_lt_thresh", "gt_frac", "lt_frac")
    return(stats_df)
}

load_xpnsl_norm <- function(path){
    chrom_name <- str_extract(basename(path), pattern = '.+(?=_Urban)')
    perm <- str_extract(basename(path), pattern =  '(?<=Rural_)\\d+(?=_permuted)')
    df <- suppressMessages(read_delim(path, delim = "\t")) %>%
        mutate(Chr = chrom_name, iter = perm)
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
    group_split(Chr, iter) %>%
    purrr::map_dfr(., calculate_windowed_stats, window_size = window_size, step = step, thresh = thresh)

write_delim(xpnsl_windows, snakemake@output[["xpnsl_df"]], delim = "\t")
