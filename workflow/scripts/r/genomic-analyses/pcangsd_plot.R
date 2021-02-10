# Script to perform PCA on covariance in allele frequencies at 4-fold sites

# Load packages
library(tidyverse)
library(wesanderson)

# Load covariance matrix
cov_mat <- read_delim('data/pcangsd/allSamples_4fold_maf0.05_pcangsd.cov', 
                      col_names = FALSE, delim = ' ') %>% 
  as.matrix()

# Load sample order
habitat_info <- read_delim('resources/sequencedPlants_phenotypesHabitat.txt', 
                           delim = '\t') %>% 
  dplyr::select(Sample, Habitat, Population, Plant)
samples <- read_table('data/angsd_sample_order.txt', col_names = FALSE) %>% 
  rename('Sample' = 'X1') %>%
  left_join(., habitat_info, by = 'Sample')
  
# Extract eigenvectors and create df
summary(princomp(cov_mat))
eigenvectors <- eigen(cov_mat)
eigen_df <- eigenvectors$vectors %>% 
  as.data.frame() %>% 
  dplyr::select(V1, V2, V3, V4) %>% 
  rename('PC1' = 'V1',
         'PC2' = 'V2', 
         'PC3' = 'V3',
         'PC4' = 'V4') %>% 
  bind_cols(., samples)

# Plot
cols <- wes_palette("Darjeeling1", n = 3, type = 'discrete')
ggplot(eigen_df, aes(x = PC1, y = PC2, color = Habitat)) +
  geom_point(size = 2.5) + 
  scale_color_manual(values = cols) + 
  theme_classic() + 
  xlab('PC1 (4.0%)') + ylab('PC2 (2.6%)') +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))

cols <- wes_palette("Darjeeling1", n = 5, type = 'discrete')
cols <-c(cols[1], cols[3])
eigen_df %>% 
  filter(PC1 < -0.2 | PC1 > 0.2 | PC2 > 0.19 | Sample %in% c('s_42_1', 's_42_20')) %>% 
  ggplot(., aes(x = PC1, y = PC2, color = Habitat)) +
  geom_point(size = 2.5) + 
  scale_color_manual(values = cols) + 
  theme_classic() + 
  xlab('PC1 (4.0%)') + ylab('PC2 (2.6%)') +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))

# Plot
cols <- wes_palette("Darjeeling1", n = 3, type = 'discrete')
ggplot(eigen_df, aes(x = PC3, y = PC4, color = Habitat)) +
  geom_point(size = 2.5) + 
  scale_color_manual(values = cols) + 
  theme_classic() + 
  xlab('PC3 (2.3%)') + ylab('PC4 (2.0%)') +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))

eigen_df %>% 
  filter(PC3 < -0.25 | PC3 > 0.25 | PC4 > 0.25 | PC4 < -0.25)

#####################################
#### EXCLUDING WEIRD INDIVIDUALS ####
#####################################

# Load covariance matrix
cov_mat_red <- read_delim('~/singularity-dev/outlierTest_pcangsd.cov', 
                      col_names = FALSE, delim = ' ') %>% 
  as.matrix()

samples_red <- samples %>% 
  filter(!(Sample %in% c('s_97_14', 's_37_18', 's_42_9', 's_42_20', 's_83_11', 's_37_11')))

# Extract eigenvectors and create df
summary(princomp(cov_mat_red))
eigenvectors_reduced <- eigen(cov_mat_red)
eigen_df_reduced <- eigenvectors_reduced$vectors %>% 
  as.data.frame() %>% 
  dplyr::select(V1, V2) %>% 
  rename('PC1' = 'V1',
         'PC2' = 'V2') %>% 
  bind_cols(., samples_red)

# Plot
cols <- wes_palette("Darjeeling1", n = 38, type = 'continuous')
ggplot(eigen_df_reduced, aes(x = PC1, y = PC2, color = factor(Population))) +
  geom_point(size = 2.5) + 
  scale_color_manual(values=cols) +
  theme_classic() + 
  xlab('PC1 (3.3%)') + ylab('PC2 (2.2%)') +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))

eigen_df_reduced %>% 
  distinct(Population)
