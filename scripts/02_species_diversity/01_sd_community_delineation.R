

## load ----
data_PA <- readRDS("data/PA_Mat_GaspObis.RDS")

list_communities <- list()


## subset PA to family ----

# subset to family of interest
fam = "Lutjanidae"

# keep only species in target family, and for which we have presence data
list_species <- 
  data_species %>% 
  filter(family %in% fam) %>%
  filter(species %in% colnames(data_PA)) %>% 
  pull(species) # keep only species info

# create PA matrix for those species
list_communities[[fam]] <- list_species
  # data_PA %>% 
  # select(all_of(c("Longitude", "Latitude", list_species)))





## subset PA by depth range ----

# create cateogries of depth
depth_list <- c(0, 30, 150, 300, 1000, 4000, 9000)
depth_categories <- c("Shallow", "Mesophotic", "Rariphotic", "Mesopelagic", "Bathypelagic") 
  
data_depth <- data.frame(matrix(NA, ncol = length(depth_list), nrow = nrow(data_species)))
rownames(data_depth) <- data_species$species
colnames(data_depth) <- depth_list

for (sp in 1: length(data_species$species)){
  species <- data_species[sp, "species"]
  indice_min <- findintervals(data_species[sp, "depth_min"], depth_list)
  indice_max <- findintervals(data_species[sp, "depth_max"], depth_list)
  data_depth[species, indice_min:indice_max] <- 1
}

data_depth[is.na(data_depth)] <- 0
data_depth <- data_depth[, 1:5]
colnames(data_depth) <- depth_categories

write.csv(data_depth, "intermediate/02_species_diversity/species_depth_categories.csv", row.names = F, quote = F)




for (depth in depth_categories){
  
  #subset species from range
  list_species <-
    data_depth %>% 
    rownames_to_column(var = "species") %>%
    filter(.data[[depth]] == 1) %>% # select species present in the depth category
    filter(species %in% colnames(data_PA)) %>% 
    pull(species) # keep species name
  
  list_communities[[depth]] <- list_species

}






## export ----
saveRDS(list_communities, "intermediate/02_species_diversity/List_community.RDS")



