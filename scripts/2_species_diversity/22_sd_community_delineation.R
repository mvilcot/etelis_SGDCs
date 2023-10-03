
# ---- load presence data ----
data_PA <- readRDS("data/PA_Mat_GaspObis.RDS")


# ---- community taxo ----
## get taxo group names ----
subfam <-
  data_taxo %>% 
  filter(Species == "Etelis_coruscans") %>% 
  pull(Subfamily)

fam <-
  data_taxo %>% 
  filter(Species == "Etelis_coruscans") %>% 
  pull(Family)

ord  <-
  data_taxo %>% 
  filter(Species == "Etelis_coruscans") %>% 
  pull(Order)

cla  <-
  data_taxo %>% 
  filter(Species == "Etelis_coruscans") %>% 
  pull(Class)



## keep only species in target group, and for which we have presence data ----
list_communities <- list()

# subfamily
list_species <- 
  data_taxo %>% 
  filter(Subfamily %in% subfam) %>%
  filter(Species %in% colnames(data_PA)) %>% 
  pull(Species) # keep only species info

list_communities[[subfam]] <- list_species

# family
list_species <- 
  data_taxo %>% 
  filter(Family %in% fam) %>%
  filter(Species %in% colnames(data_PA)) %>% 
  pull(Species) # keep only species info

list_communities[[fam]] <- list_species

# order
list_species <- 
  data_taxo %>% 
  filter(Order %in% ord) %>%
  filter(Species %in% colnames(data_PA)) %>% 
  pull(Species) # keep only species info

list_communities[[ord]] <- list_species

# class
list_species <- 
  data_taxo %>% 
  filter(Class %in% cla) %>%
  filter(Species %in% colnames(data_PA)) %>% 
  pull(Species) # keep only species info

list_communities[[cla]] <- list_species



## export ----
saveRDS(list_communities, "intermediate/2_species_diversity/List_community_taxonomy.RDS")




# ---- community depth X taxo ----
comm_delin = "taxonomy"
list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

depth_Ec_min <- 45
depth_Ec_max <- 400


## 1.crosses range ----
# select species whose depth range intersects that of Etelis coruscans
data_depth$crosses45_400 <- 
  ifelse((data_depth$depth_max <= depth_Ec_min) |  # if depth_max lower than shallower
           (data_depth$depth_min >= depth_Ec_max), # or depth_min higher than deeper
         "no",                                     # than not crossing Etelis coruscans range
         "yes")

data_depth <- 
  data_depth %>% 
  relocate(crosses45_400, .after = depth_max) 

list_crosses45_400 <-
  data_depth %>% 
  filter(crosses45_400 == "yes") %>% 
  pull(species)


#subset species that crosses 45-400m range
list_communitiesDEPTH <- list()
for(comm in names(list_communities)){
  list_communitiesDEPTH[[comm]] <- 
    list_communities[[comm]][list_communities[[comm]] %in% list_crosses45_400]
}

# export
saveRDS(list_communitiesDEPTH, "intermediate/2_species_diversity/List_community_taxonomy_depth1_crosses45-400m.RDS")



## 2.contains range ----
# select species whose depth range contains (or is within) that of Etelis coruscans
data_depth$contains45_400 <-
  ifelse(((data_depth$depth_max <= depth_Ec_max) &  
            (data_depth$depth_min >= depth_Ec_min))| # if range is within Etelis corucsans'
           ((data_depth$depth_max >= depth_Ec_max) &  
              (data_depth$depth_min <= depth_Ec_min)),  # OR if overlaps Etelis corucsans" range
         "yes",                                      # than keep species
         "no")

data_depth <- 
  data_depth %>% 
  relocate(contains45_400, .after = depth_max) 

list_contains45_400 <-
  data_depth %>% 
  filter(contains45_400 == "yes") %>% 
  pull(species)


#subset species within or containing 45-400m range
list_communitiesDEPTH <- list()
for(comm in names(list_communities)){
  list_communitiesDEPTH[[comm]] <- 
    list_communities[[comm]][list_communities[[comm]] %in% list_contains45_400]
}

# export
saveRDS(list_communitiesDEPTH, "intermediate/2_species_diversity/List_community_taxonomy_depth2_contains45-400m.RDS")



## 3.within range ----
# select species whose depth range is contains within that of Etelis coruscans
data_depth$within45_400 <-
  ifelse(((data_depth$depth_max <= depth_Ec_max) &  
            (data_depth$depth_min >= depth_Ec_min)), # if range is within Etelis corucsans'
         "yes",                                      
         "no")

data_depth <- 
  data_depth %>% 
  relocate(within45_400, .after = depth_max) 

list_within45_400 <-
  data_depth %>% 
  filter(within45_400 == "yes") %>% 
  pull(species)


#subset species within 45-400m range
list_communitiesDEPTH <- list()
for(comm in names(list_communities)){
  list_communitiesDEPTH[[comm]] <- 
    list_communities[[comm]][list_communities[[comm]] %in% list_within45_400]
}

# export
saveRDS(list_communitiesDEPTH, "intermediate/2_species_diversity/List_community_taxonomy_depth3_within45-400m.RDS")




# ---- community envt X taxo ----
data_depth$reef_associated <-
  ifelse(data_depth$env_2 == "reef-associated",
         "yes",                                      
         "no")

data_depth <- 
  data_depth %>% 
  relocate(reef_associated, .after = depth_max) 

list_reef_associated <-
  data_depth %>% 
  filter(reef_associated == "yes") %>% 
  pull(species)


#subset species within 45-400m range
list_communitiesENV <- list()
for(comm in names(list_communities)){
  list_communitiesENV[[comm]] <- 
    list_communities[[comm]][list_communities[[comm]] %in% list_reef_associated]
}

# export
saveRDS(list_communitiesENV, "intermediate/2_species_diversity/List_community_taxonomy_env_reef-associated.RDS")




# ---- *** DRAFTS ----
## ---- *** subset PA by phylogenetic scale ----
# 
# # compute distance from Etelis_coruscans
# tree_all <- fishtree_phylogeny()
# dist_phylo <- cophenetic.phylo(tree_all)
# dist_Ecoruscans <- 
#   data.frame(phylodist = dist_phylo[, "Etelis_coruscans"]) %>% 
#   rownames_to_column("species")
# 
# write.csv(dist_Ecoruscans, "intermediate/2_species_diversity/Phylogenetic_distance_to_Etelis_coruscans.csv",
#           row.names = F, quote = F)
# 
# # keep only species in target community, and for which we have presence data
# list_communities <- list()
# dist_chronogram <- c(100,120,140,160,200,300,400,1000)
# # dist_phylogram <- seq(from = 0.2, to = 4, by = 0.1)
# 
# for (d in dist_chronogram){
#   
#   comm <- paste0("phylodist_", d)
#   print(comm)
#   
#   list_communities[[comm]] <-
#     dist_Ecoruscans %>% 
#     filter(phylodist < d) %>% 
#     filter(species %in% colnames(data_PA)) %>% 
#     pull(species)
# }
# 
# saveRDS(list_communities, "intermediate/2_species_diversity/List_community_phylogenetic_distance.RDS")
# 
# 
# 
## ---- *** subset PA by taxonomic scales from Zurich data ----
# list_communities <- list()
# 
# # subset to family of interest
# fam = "Lutjanidae"
# ord = "Lutjaniformes"
# cla = "Actinopteri"
# 
# # keep only species in target family, and for which we have presence data
# list_species <- 
#   data_depth2 %>% 
#   filter(family %in% fam) %>%
#   filter(species %in% colnames(data_PA)) %>% 
#   pull(species) # keep only species info
# 
# list_communities[[fam]] <- list_species
# 
# # keep only species in target family, and for which we have presence data
# list_species <- 
#   data_depth2 %>% 
#   filter(order %in% ord) %>%
#   filter(species %in% colnames(data_PA)) %>% 
#   pull(species) # keep only species info
# 
# list_communities[[ord]] <- list_species
# 
# # keep only species in target family, and for which we have presence data
# list_species <- 
#   data_depth2 %>% 
#   filter(class %in% cla) %>%
#   filter(species %in% colnames(data_PA)) %>% 
#   pull(species) # keep only species info
# 
# list_communities[[cla]] <- list_species
# 
# # save
# saveRDS(list_communities, "intermediate/2_species_diversity/List_community_taxonomic_scale_datasp2.RDS")
# 
#
#
#
#
## ---- *** subset PA by depth range ----
# 
# # -- subset PA to family --
# list_communities <- list()
# fam = "Lutjanidae"
# 
# # keep only species in target family, and for which we have presence data
# list_species <- 
#   data_depth %>% 
#   filter(family %in% fam) %>%
#   filter(species %in% colnames(data_PA)) %>% 
#   pull(species) # keep only species info
# 
# list_communities[[fam]] <- list_species
# 
# 
# # -- create cateogries of depth --
# depth_list <- c(0, 30, 150, 300, 1000, 4000, 9000)
# depth_categories <- c("Shallow", "Mesophotic", "Rariphotic", "Mesopelagic", "Bathypelagic") 
#   
# data_depth <- data.frame(matrix(NA, ncol = length(depth_list), nrow = nrow(data_depth)))
# rownames(data_depth) <- data_depth$species
# colnames(data_depth) <- depth_list
# 
# for (sp in 1: length(data_depth$species)){
#   species <- data_depth[sp, "species"]
#   indice_min <- findintervals(data_depth[sp, "depth_min"], depth_list)
#   indice_max <- findintervals(data_depth[sp, "depth_max"], depth_list)
#   data_depth[species, indice_min:indice_max] <- 1
# }
# 
# data_depth[is.na(data_depth)] <- 0
# data_depth <- data_depth[, 1:5]
# colnames(data_depth) <- depth_categories
# 
# # write.csv(data_depth, "intermediate/2_species_diversity/species_depth_categories.csv", 
# #           row.names = T, quote = F)
# 




