# Load required packages
library(tidyverse)

# Load all data
arg_fst_df <- suppressMessages(snakemake@input[["arg_fst"]] %>% purrr::map_dfr(., read_csv))

gt_fst_df <- suppressMessages(read_delim(snakemake@input[["gt_fsts"]], delim='\t') %>% filter(pop1 == "Urban" & pop2 == "Rural"))

cols <- c("chrom", "pos", "fst_num", "fst_denom")
sfs_fst_df <- suppressMessages(read_delim(snakemake@input[["sfs_fsts"]], delim = "\t", col_names = cols))

load_region_files <- function(path){
    region_id = str_split(basename(path), "\\.", simplify=T)[1,4]
    dat <- suppressMessages(
        read_delim(path, delim="\t", col_names=c("chrom", "start", "end")) %>% 
        mutate(region_id = region_id)
    )
    return(dat)
}

regions_df <- snakemake@input[["regions"]] %>%
    purrr::map_df(., load_region_files)

merge_gt_fst_vals <- function(df){

    region <- df %>% pull(regionID) %>% unique()
    region_df <- regions_df %>% filter(region_id == region) 
    start <- region_df$start + 1
    end <- region_df$end
    window_size <- snakemake@params[["window_size"]]
    
    arg_fsts <- df %>% 
        dplyr::select(regionID, arg_win_start, arg_win_end, arg_branch_fst, arg_site_fst) %>% 
        rename("window_pos_1" = "arg_win_start",
               "window_pos_2" = "arg_win_end")
    gt_fsts <- gt_fst_df %>% 
        filter(window_pos_1 >= start & window_pos_2 <= end)

    all_fsts <- gt_fsts %>% 
        left_join(arg_fsts, by = c("window_pos_1", "window_pos_2"))
    return(all_fsts)
}

arg_gt_fsts_df <- arg_fst_df %>%
    group_split(regionID) %>% 
    purrr::map_dfr(., merge_gt_fst_vals)

merge_sfs_fst_vals <- function(starts, ends){

    sfs_hudson_fst <- sfs_fst_df %>% 
        filter(pos >= starts & pos <= ends) %>% 
        summarise(sfs_hudson_fst = sum(fst_num) / sum(fst_denom)) %>% 
        pull(sfs_hudson_fst)
    return(sfs_hudson_fst)
}

all_fsts_df <- arg_gt_fsts_df %>% 
    mutate(sfs_hudson_fst = purrr::map2_dbl(.$window_pos_1, .$window_pos_2, merge_sfs_fst_vals, .progress = T))

write_delim(all_fsts_df, snakemake@output[[1]], delim="\t")

