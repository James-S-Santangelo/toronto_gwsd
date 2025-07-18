# Scripts to generate files used to create the map of sampled populations with HII by habitat inset
# The map itself will be created in QGIS

library(tidyverse)

plants <- read_delim("resources/sequencedPlants_phenotypesHabitat.txt", delim = '\t') %>% 
  dplyr::select(Habitat, Population, Plant)
all_plants <- read_csv("resources/reference/allPlants_Toronto.csv") %>% 
  dplyr::select(Population, Plant, Transect, Lat.pop, Long.pop) %>% 
  mutate(Transect = case_when(Transect == "A" ~ "North",
                              Transect == "B" ~ "West",
                              Transect == "C" ~ "East"))
transects <- all_plants %>% 
  dplyr::select(-Plant) %>% 
  distinct()
hii <- read_csv("resources/reference/Toronto_hii.csv") %>% 
  dplyr::select(population, hii) %>% 
  rename("Population" = "population") %>% 
  left_join(transects, by = "Population")



# Create dataframe for plotting mean HII by Habitat. This is inset to map
inset_df <- plants %>% 
  dplyr::select(-Plant) %>% 
  left_join(., hii, by = "Population") %>% 
  mutate(Transect = ifelse(Habitat == "Urban", "All", Transect)) %>% 
  distinct() %>% 
  group_by(Habitat, Transect) %>% 
  summarise(mean = mean(hii),
            sd = sd(hii),
            se = sd / sqrt(n()))
inset_df

# Create inset plot
cols_hab <- c("#007243", "#914205", "#003876")
dodge <- position_dodge(0.75)
err_width <- ifelse(Habitat == "Urban", 0.1, 0.3)
map_inset <- ggplot(data = inset_df, aes(x = Habitat, y = mean, group = Transect)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = Habitat), size = 1,
                width = ifelse(inset_df$Habitat == "Urban", 0.1, 0.3), show.legend = F, position=dodge) +
  geom_point(aes(shape = Habitat, color = Habitat), size = 7, show.legend = F,
             position=dodge) +
  scale_color_manual(values = cols_hab) +
  theme_classic() +
  ylab("Mean human influence index (HII)") +
  geom_text(aes(label = Transect, y =  mean + se), vjust = -0.75, 
            position=dodge) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.title.x = element_blank())
map_inset

figures_dir <- "../figures_tables/map_with_inset"
dir.create(figures_dir, showWarnings = T)
ggsave(paste0(figures_dir, "/map_inset.pdf"), plot = map_inset, device = "pdf", width = 8,
       height = 8, units = "in", dpi = 600)

# Create dataframe that will be used to create map
map_df <- plants %>% 
  left_join(., all_plants, by = c("Population", "Plant")) %>%
  left_join(., hii, by = "Population") %>% 
  group_by(Population) %>% 
  mutate(n = n())
map_df

write_csv(map_df, paste0(figures_dir, "/map_df.csv"))
