# Script to write windowed fst/theta and xpnsl dataframes

# Load required packages
library(tidyverse)

################
#### XP-nSL ####
################

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
    perm <- df %>% pull(iter) %>% unique()
    winStarts <- seq(from = min(df$pos), to = max(df$pos), by = step)
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
    purrr::map_dfr(., load_xpnsl_norm, .progress=TRUE) %>%
    rename("normxpnsl" = "normxpehh") %>%
    dplyr::select(-id)

# Function to assign outliers to XP-nSL windows
assign_outliers <- function(df){

    nSites_thresh <- as.numeric(snakemake@params[['nSites_xpnsl']]) # Require at least this many site in a window
    df_filt <- df %>%
        mutate_at(vars(-("Chr")), as.numeric) %>% 
        filter(n >= nSites_thresh)

    # Get critical values for mean XP-nSL score and proportions greater or lesser than 2 and -2, respectively
    xpnsl_score_quant_filt <- quantile(df_filt %>% pull(mean), probs = c(0.01, 0.99))
    xpnsl_gtprop_quant_filt <- quantile(df_filt %>% pull(gt_frac), probs = 0.99)
    xpnsl_ltprop_quant_filt <- quantile(df_filt %>% pull(lt_frac), probs = 0.99)

    # Identify outliers and add as categorical variable to windows dataframe
    df_filt <- df_filt %>%
        mutate(xpnsl_score_outlier = ifelse(mean <= xpnsl_score_quant_filt[1] | mean >= xpnsl_score_quant_filt[2], 1, 0),
               xpnsl_gtprop_outlier = ifelse(gt_frac >= xpnsl_gtprop_quant_filt, 1, 0),
               xpnsl_ltprop_outlier = ifelse(lt_frac >= xpnsl_ltprop_quant_filt, 1, 0),
               direction = case_when(xpnsl_score_outlier == 1 & mean > 0 & xpnsl_gtprop_outlier == 1 ~ 'Urban sel',
                                     xpnsl_score_outlier == 1 & mean < 0 & xpnsl_ltprop_outlier == 1 ~ 'Rural sel',
                                     TRUE ~ 'Not outlier')) %>% 
    mutate(prop_outlier = case_when(direction == 'Urban sel' ~ gt_frac,
                                    direction == 'Rural sel' ~ lt_frac,
                                    TRUE ~ NA))

    return(df_filt)

}

# Calculate windowed xpnsl
window_size <- as.numeric(snakemake@params["winsize"])
step <- as.numeric(snakemake@params["winsize"])
thresh <- 2
xpnsl_windows <- xpnsl_norm_df %>%
    group_split(Chr, iter) %>%
    purrr::map(., calculate_windowed_stats, window_size = window_size, step = step, thresh = thresh, .progress=TRUE) %>%
    purrr::map_dfr(., assign_outliers, .progress=TRUE)

write_delim(xpnsl_windows, snakemake@output[["xpnsl_df"]], delim = "\t")
