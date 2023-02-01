## ---- genetic data ----
data_samples <- read.csv("data/metadata_samples.csv")
data_sites <- read.csv("data/metadata_sites.csv")

data_samples <- 
  data_samples %>% 
  left_join(data_sites, by = c("station", "site"))


## ---- species data ----
data_Etelis <- readRDS("data/Presence_data_Fishbase_Etelis_coruscans.RDS")
data_species <- read.csv("data/data_species_depth_range_teleo.csv")
# data_species <- read.csv("data/data_species.csv")





# # check species database
# temp1 <- 
#   data_depth %>% 
#   dplyr::filter(family == "Lutjanidae")
# 
# temp2 <- 
#   data_species %>% 
#   dplyr::filter(family == "Lutjanidae")
# 
# species1 <- temp1$species
# species2 <- temp2$species
# 
# species1[!(species1 %in% species2)]
# species2[!(species2 %in% species1)]

