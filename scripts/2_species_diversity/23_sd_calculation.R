
## >>>>> !!!!! CHANGE SCRIPT AND MAKE IT AS A FUNCTION, depending of community delineation ####

# This script calculates species richness per site and jaccard indices between sites


# ---- load data ----
## presence data
PAstation <- readRDS("intermediate/2_species_diversity/PA_Mat_GaspObis_station.RDS")
PAsite <- readRDS("intermediate/2_species_diversity/PA_Mat_GaspObis_site.RDS")
PAall <- readRDS("intermediate/2_species_diversity/PA_Mat_GaspObis_allstations.RDS")

## communtity delineation
# comm_delin = "taxonomy"
# comm_delin = "taxonomy_depth1_crosses45-400m"
# comm_delin = "taxonomy_depth2_contains45-400m"
# comm_delin = "taxonomy_depth3_within45-400m"
comm_delin = "taxonomy_env_reef-associated"

list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))
names(list_communities)

# ---- initiate ----
list_sd_gamma <- list()
list_sd_alpha_site <- list()
list_sd_alpha_station <- list()
list_sd_beta_site_pair <- list()
list_sd_beta_site_multi <- list()
list_sd_beta_station_pair <- list()
list_sd_beta_station_multi <- list()


for (comm in names(list_communities)){
  
  print(comm)

  # ---- subset to community ----
  PAstation_comm <-
    PAstation %>% 
    dplyr::select(all_of(c("station", list_communities[[comm]])))
  
  PAsite_comm <-
    PAsite %>% 
    dplyr::select(all_of(c("site", list_communities[[comm]])))
  
  PAall_comm <-
    PAall %>% 
    dplyr::select(all_of(list_communities[[comm]]))
  
  # ---- alpha diversity ----
  # Species richness at each station
  list_sd_alpha_station[[comm]] <-
    PAstation_comm %>%
    pivot_longer(cols = -c(station), names_to = "species") %>%
    group_by(station) %>% 
    summarise(community = comm, 
              richness_station = sum(value), .groups = "keep")
  
  # Species richness at each site
  list_sd_alpha_site[[comm]] <-
    PAsite_comm %>%
    pivot_longer(cols = -c(site), names_to = "species") %>%
    group_by(site) %>% 
    summarise(community = comm, 
              richness_site = sum(value), .groups = "keep")
  
  # ---- gamma diversity ----
  list_sd_gamma[[comm]] <-
    data.frame(community = comm) %>% 
    cbind(richness_allstations = rowSums(PAall_comm))
    

  # ---- beta diversity ----

  # beta diversity station level
  PAstation_comm <- 
    column_to_rownames(PAstation_comm, var = "station")
  
  list_sd_beta_station_multi[[comm]] <- 
    beta.multi(PAstation_comm, index.family="jaccard")
  
  list_sd_beta_station_pair[[comm]] <- 
    beta.pair(PAstation_comm, index.family="jaccard")
  
  # beta diversity site level
  PAsite_comm <- 
    column_to_rownames(PAsite_comm, var = "site")
  
  list_sd_beta_site_multi[[comm]] <- 
    beta.multi(PAsite_comm, index.family="jaccard")
  
  list_sd_beta_site_pair[[comm]] <- 
    beta.pair(PAsite_comm, index.family="jaccard")

}





# ---- collate ----
# alpha station
sd_alpha_station <- 
  do.call(rbind.data.frame, list_sd_alpha_station)

# alpha site
sd_alpha_site <- 
  do.call(rbind.data.frame, list_sd_alpha_site)

# alpha mean station
sd_alpha_mean_station <- 
  sd_alpha_station %>%
  group_by(community) %>% 
  summarise(richness_station_mean = mean(richness_station), 
            .groups = "keep")

# alpha mean site
sd_alpha_mean_site <- 
  sd_alpha_site %>%
  group_by(community) %>% 
  summarise(richness_site_mean = mean(richness_site), 
            .groups = "keep")

# gamma
sd_gamma <- 
  do.call(rbind.data.frame, list_sd_gamma) %>% 
  as_tibble()


# beta
sd_beta_station <- do.call(rbind.data.frame, list_sd_beta_station_multi) 
colnames(sd_beta_station) <- 
  paste0(colnames(sd_beta_station), "_station")
sd_beta_station <-
  sd_beta_station %>% 
  rownames_to_column(var = "community")

sd_beta_site <- do.call(rbind.data.frame, list_sd_beta_site_multi) 
colnames(sd_beta_site) <- 
  paste0(colnames(sd_beta_site), "_site")
sd_beta_site <-
  sd_beta_site %>% 
  rownames_to_column(var = "community")



# ---- merge ----
sd_global <-
  sd_gamma %>% 
  merge(sd_alpha_mean_site) %>% 
  merge(sd_alpha_mean_station) %>% 
  merge(sd_beta_site) %>% 
  merge(sd_beta_station)



# ---- export ----
sd_global %>% write.csv(paste0("results/2_species_diversity/sd_table_global_", comm_delin, ".csv"), row.names = FALSE)
sd_alpha_site %>%  write.csv(paste0("results/2_species_diversity/sd_table_site_", comm_delin, ".csv"), row.names = FALSE)
sd_alpha_station %>% write.csv(paste0("results/2_species_diversity/sd_table_station_", comm_delin, ".csv"), row.names = FALSE)
list_sd_beta_site_pair %>% saveRDS(paste0("results/2_species_diversity/sd_list_pairwise_site_", comm_delin, ".RDS"))
list_sd_beta_station_pair %>% saveRDS(paste0("results/2_species_diversity/sd_list_pairwise_station_", comm_delin, ".RDS"))







# ---- *** DRAFTS ----
## ---- *** Compare Lutjanidae dataset ----
# 
# lutjanidae_list <- list()
# comm = "Lutjanidae"
# for (comm_delin in c("depth_category", "taxonomic_scale_datasp2", "taxonomic_scale_Fishbase")){
#   
#   list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))
#   
#   ## subset to community
#   lutjanidae_list[[comm_delin]] <-
#     list_communities[[comm]]
#   
# }
# 
# tree_lutj  <- fishtree_phylogeny(rank = "Lutjanidae")
# lutjanidae_list[["Rabosky"]] <- tree_lutj$tip.label
# 
# 
# library(qdapTools)
# lutjanidae_df <- t(as.matrix(mtabulate(lutjanidae_list)))
# write.csv(lutjanidae_df, "intermediate/2_species_diversity/Lutjanidae_datasets_comparison.csv", quote = F, row.names = T)
# 
# 
# 
# 
# 
# 
# 
## ---- *** Compare Lutjanidae dataset in PA ----
# 
# lutjanidae_PA <- list()
# comm = "Lutjanidae"
# for (comm_delin in c("depth_category", "taxonomic_scale_datasp2", "taxonomic_scale_Fishbase")){
#   
#   list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))
#   
#   ## subset to community
#   PAall_comm <-
#     PAall %>% 
#     dplyr::select(all_of(list_communities[[comm]]))
#   
#   ## gamma diversity
#   data.frame(community = comm) %>% 
#     cbind(richness_allstations = rowSums(PAall_comm))
#   
#   ## compare Lutjanidae database
#   if (comm == "Lutjanidae"){
#     lutjanidae_PA[[comm_delin]] <- colnames(PAall_comm[, colSums(PAall_comm) != 0])
#   }
# }
# 
# 
# 
# tree_lutj  <- fishtree_phylogeny(rank = "Lutjanidae")
# 
# # temp <- data.frame(colnames(PAall_comm[, colSums(PAall_comm) != 0]))
# # lutjanidae_df <- as.data.frame(do.call(cbind, lutjanidae_PA))
# # write.csv(lutjanidae_df, "Lutjanidae_presence_all_sites_datasets_comparison.csv", quote = F, row.names = F)
# 
# library(qdapTools)
# lutjanidae_df <- t(as.matrix(mtabulate(lutjanidae_PA)))
# write.csv(lutjanidae_df, "intermediate/2_species_diversity/Lutjanidae_datasets_comparison_presence_all_sites_PA.csv", quote = F, row.names = T)
# 
