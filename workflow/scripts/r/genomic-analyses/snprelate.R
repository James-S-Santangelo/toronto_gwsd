# Script to look at pairwise relatedness within populations

# Load packages
library(tidyverse)

###################
#### FUNCTIONS ####
###################

load_ngsrelate_data <- function(path){
  
  # Get name of folder with parameter combinations
  dirname <- dirname(path)
  
  # Read in NGSrelate data
  full_path <- paste0(inpath, '/', path)
  dat <- read_delim(full_path, delim= '\t', col_names = TRUE) %>%
    as.data.frame() %>% 
    mutate(pop = dirname)
  return(dat)
}

##################
#### ANALYSIS ####
##################

inpath <- 'data/ngsrelate/'
ngsrelate_df <- list.files(inpath, recursive = TRUE) %>% 
  map_dfr(., load_ngsrelate_data)
  
ngsrelate_df %>% 
  dplyr::select(a, b, nSites, rab, Fa, Fb, pop) %>% 
  filter(pop == 'pop7')
  arrange(desc(rab)) %>% 
  View()

