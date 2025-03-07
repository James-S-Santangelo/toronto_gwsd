# Load required packages
library(tidyverse)

calculate_windowed_ld <- function(df, window_size, step, thresh){
    chrom <- df %>% pull(CHR_A) %>% unique()
    pos_vector <- df %>% pull(BP_A) %>% unique()
    winStarts <- seq(from = 0, to = max(pos_vector) + window_size, by = step)
    mat <- matrix(0, nrow = length(winStarts), ncol = 9)
    for(i in 1:length(winStarts)){
        start <- winStarts[i]
        end <- start + step
        df_filt <- df %>% filter(BP_A >= start & BP_A < end)
        winID <- i
        winCenter <- start + (step / 2)
        mean <- suppressWarnings(mean(df_filt$R2))
        max <- suppressWarnings(max(df_filt$R2))
        min <- suppressWarnings(min(df_filt$R2))
        n <- length(pos_vector)
        stats <- c(chrom, winID, start, end, winCenter, mean, max, min, n)
        mat[i, ] <- stats
    }
    stats_df <- as.data.frame(mat)
    names(stats_df) <- c("Chr", "winID", "start", "end", "winCenter", "mean", "max", "min", "n")
    return(stats_df)
}


habitat <- snakemake@wildcards[["habitat"]]
window_size <- as.numeric(snakemake@params[["winsize"]])
step <- as.numeric(snakemake@params[["winsize"]])

# Load Pairwise LD data
ld_df <- read_table(snakemake@input[["ld"]]) %>% dplyr::select(-X8)

# Calculate windowed haplotype stats
ld_windowed <- calculate_windowed_ld(ld_df,
                                     window_size = window_size,
                                     step = step,
                                     thresh = thresh) %>%
    mutate(habitat = habitat)

write_delim(ld_windowed, snakemake@output[["win_ld"]], delim = "\t")
