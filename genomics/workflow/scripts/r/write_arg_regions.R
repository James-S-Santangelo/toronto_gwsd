library(tidyverse)

# Load dataframes with windowed sweep statistics
win_xpnsl_df <- read_delim(snakemake@input[['xpnsl']], delim = '\t')

urban_rural_windows <- win_xpnsl_df %>% 
    filter(direction %in% c('Urban sel', 'Rural sel'))
unselected_window <- win_xpnsl_df %>% 
    filter(direction == "Not outlier") %>% 
    sample_n(120)

all_regions <- bind_rows(urban_rural_windows, unselected_window) %>% 
    mutate(region = paste0(Chr, ":", start, "-", end),
           regionID = 1:n()) %>% 
    dplyr::select(Chr, start, end, winID, region, regionID, direction)

write_delim(x = all_regions, file = snakemake@output[["arg_regions"]], delim = "\t")
