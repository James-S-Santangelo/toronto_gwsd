# Script to write windowed fst/theta and xpnsl dataframes

# Load required packages
library(tidyverse)

################
#### XP-nSL ####
################

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

calculate_windowed_hapstat <- function(df, var, window_size, step, thresh){
    chrom <- df %>% pull(Chr) %>% unique()
    habitat <- df %>% pull(habitat) %>% unique()
    winStarts <- seq(from = 0, to = max(df$pos) + window_size, by = step)
    mat <- matrix(0, nrow = length(winStarts), ncol = 11)
    for(i in 1:length(winStarts)){
        start <- winStarts[i]
        end <- start + step
        df_filt <- df %>% filter(pos >= start & pos < end)
        winID <- i
        winCenter <- start + (step / 2)
        stat_vals <- df_filt %>% pull(!!enquo(var))
        mean <- mean(stat_vals)
        max <- max(stat_vals)
        n <- nrow(df_filt)
        num_gt_thresh <- sum(stat_vals > thresh)
        gt_frac <- sum(num_gt_thresh) / n
        stats <- c(chrom, habitat, winID, start, end, winCenter, mean, max, n, num_gt_thresh, gt_frac)
        mat[i, ] <- stats
    }
    stats_df <- as.data.frame(mat)
    names(stats_df) <- c("Chr", "habitat", "winID", "start", "end", "winCenter", "mean", "max", "n", "num_gt_thresh", "gt_frac")
    return(stats_df)
}


window_size <- as.numeric(snakemake@params["winsize"])
step <- as.numeric(snakemake@params["winsize"])
thresh <- 2

if(snakemake@wildcards[["stat"]] == "xpnsl"){

    load_xpnsl_norm <- function(path){
        chrom_name <- str_extract(basename(path), pattern = '.+(?=_Urban)')
        df <- suppressMessages(read_delim(path, delim = "\t")) %>%
            mutate(Chr = chrom_name)
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
} else {

    load_hapstat_norm <- function(path, stat) {
        chrom_name <- str_extract(basename(path), pattern = ".+(?=_[Urban|Rural])")
        habitat <- str_extract(basename(path), pattern = "(Urban|Rural)")

        if (stat == "ihh12") {
            df <- suppressMessages(read_delim(path, delim = "\t")) %>%
                mutate(Chr = chrom_name, habitat = habitat)
        } else {
            if (stat == "ihs") {
                cols <- c("id", "pos", "1_freq", "ihh1", "ihh0", "iHs", "normihs", "crit")
                df <- suppressMessages(read_delim(path, delim = "\t", col_names = cols)) %>%
                    mutate(Chr = chrom_name, habitat = habitat) %>%
                    mutate(abs_normihs = abs(normihs))
            } else if (stat == "nsl") {
                cols <- c("id", "pos", "1_freq", "sl1", "sl0", "nSL", "normnsl", "crit")
                df <- suppressMessages(read_delim(path, delim = "\t", col_names = cols)) %>%
                    mutate(Chr = chrom_name, habitat = habitat) %>%
                    mutate(abs_normnsl = abs(normnsl))
            }
        }
        return(df)

    }

    # Load raw ihh12 data
    stats_norm_df <- snakemake@input[["norm"]] %>%
        purrr::map_dfr(., load_hapstat_norm, stat=snakemake@wildcards[["stat"]]) %>%
        dplyr::select(-id)

    if (snakemake@wildcards[["stat"]] == "ihh12") {
        var <- "normihh12"
    } else if (snakemake@wildcards[["stat"]] == "ihs") {
        var <- "abs_normihs"
    } else {
        var <- "abs_normnsl"
    }
    # Calculate windowed haplotype stats
    stats_windowed <- stats_norm_df %>%
        group_split(habitat, Chr) %>%
        purrr::map_dfr(., calculate_windowed_hapstat, var=var, window_size = window_size, step = step, thresh = thresh)

    write_delim(stats_windowed, snakemake@output[["hapstats_df"]], delim = "\t")
}


