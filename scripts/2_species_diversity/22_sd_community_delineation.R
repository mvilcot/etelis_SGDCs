

## ---- load ----
data_PA <- readRDS("data/PA_Mat_GaspObis.RDS")




## ---- subset PA by phylogenetic scale ----

# compute distance from Etelis_coruscans
tree_all <- fishtree_phylogeny()
dist_phylo <- cophenetic.phylo(tree_all)
dist_Ecoruscans <- 
  data.frame(phylodist = dist_phylo[, "Etelis_coruscans"]) %>% 
  rownames_to_column("species")

write.csv(dist_Ecoruscans, "intermediate/2_species_diversity/Phylogenetic_distance_to_Etelis_coruscans.csv",
          row.names = F, quote = F)

# keep only species in target community, and for which we have presence data
list_communities <- list()
dist_chronogram <- c(100,120,140,160,200,300,400,1000)
# dist_phylogram <- seq(from = 0.2, to = 4, by = 0.1)

for (d in dist_chronogram){
  
  comm <- paste0("phylodist_", d)
  print(comm)
  
  list_communities[[comm]] <-
    dist_Ecoruscans %>% 
    filter(phylodist < d) %>% 
    filter(species %in% colnames(data_PA)) %>% 
    pull(species)
}

saveRDS(list_communities, "intermediate/2_species_diversity/List_community_phylogenetic_distance.RDS")




## ---- subset PA by taxonomic scales ----
list_communitiesFB <- list()

# subset to family of interest
subfam = "Etelinae"
fam = "Lutjanidae"
ord = "Lutjaniformes"
ord = "Eupercaria/misc"
cla = "Actinopteri"

# keep only species in target family, and for which we have presence data
list_species <- 
  data_fishbase %>% 
  filter(subfamily %in% subfam) %>%
  filter(species %in% colnames(data_PA)) %>% 
  pull(species) # keep only species info

list_communitiesFB[[subfam]] <- list_species


# keep only species in target family, and for which we have presence data
list_species <- 
  data_fishbase %>% 
  filter(family %in% fam) %>%
  filter(species %in% colnames(data_PA)) %>% 
  pull(species) # keep only species info

list_communitiesFB[[fam]] <- list_species


# keep only species in target family, and for which we have presence data
list_species <- 
  data_fishbase %>% 
  filter(order %in% ord) %>%
  filter(species %in% colnames(data_PA)) %>% 
  pull(species) # keep only species info

list_communitiesFB[[ord]] <- list_species


# keep only species in target family, and for which we have presence data
list_species <- 
  data_fishbase %>% 
  filter(class %in% cla) %>%
  filter(species %in% colnames(data_PA)) %>% 
  pull(species) # keep only species info

list_communitiesFB[[cla]] <- list_species


saveRDS(list_communitiesFB, "intermediate/2_species_diversity/List_community_taxonomic_scale_Fishbase.RDS")




## ---- subset PA by depth range ----

# -- subset PA to family --
list_communities <- list()
fam = "Lutjanidae"

# keep only species in target family, and for which we have presence data
list_species <- 
  data_species %>% 
  filter(family %in% fam) %>%
  filter(species %in% colnames(data_PA)) %>% 
  pull(species) # keep only species info

list_communities[[fam]] <- list_species


# -- create cateogries of depth --
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

# write.csv(data_depth, "intermediate/2_species_diversity/species_depth_categories.csv", 
#           row.names = T, quote = F)


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

# 
# list_communities[["submesophotic"]] <-
#   data_species %>% 
#   filter(depth_min <= 100) %>% 
#   filter(depth_max >= 350) %>% 
#   filter(species %in% colnames(data_PA)) %>% 
#   pull(species) # keep species name


# export
saveRDS(list_communities, "intermediate/2_species_diversity/List_community_depth_category.RDS")




