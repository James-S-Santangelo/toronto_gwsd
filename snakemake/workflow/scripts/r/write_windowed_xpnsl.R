# Script to write windowed xpnsl dataframes
# Run for Urban-Rural, Suburban-Rural, Urban-Suburban

# Load required packages
library(tidyverse)

calculate_windowed_xpnsl <- function(df, window_size, step, thresh){
    chrom <- df %>% pull(Chr) %>% unique()
    winStarts <- seq(from = 0, to = max(df$pos) + window_size, by = step)
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


hab_comb <- snakemake@wildcards["hab_comb"]
window_size <- as.numeric(snakemake@params["winsize"])
step <- as.numeric(snakemake@params["winsize"])
thresh <- 2

load_xpnsl_norm <- function(path){
    if (hab_comb %in% c("Urban_Rural", "Urban_Suburban")) {
        chrom_name <- str_extract(basename(path), pattern = '.+(?=_Urban)')
        df <- suppressMessages(read_delim(path, delim = "\t")) %>%
            mutate(Chr = chrom_name)
    }else if (hab_comb == "Suburban_Rural") {
        chrom_name <- str_extract(basename(path), pattern = '.+(?=_Suburban)')
        df <- suppressMessages(read_delim(path, delim = "\t")) %>%
            mutate(Chr = chrom_name)
    }
    return(df)
}

# Load raw xpnsl data
stats_norm_df <- snakemake@input[["norm"]] %>%
    purrr::map_dfr(., load_xpnsl_norm) %>%
    rename("normxpnsl" = "normxpehh") %>%
    dplyr::select(-id)

# Calculate windowed haplotype stats 
stats_windowed <- stats_norm_df %>%
    group_split(Chr) %>%
    purrr::map_dfr(., calculate_windowed_xpnsl, window_size = window_size, step = step, thresh = thresh)

write_delim(stats_windowed, snakemake@output[["hapstats_df"]], delim = "\t")
