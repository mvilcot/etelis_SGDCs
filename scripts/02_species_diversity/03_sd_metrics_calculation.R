# This script calculates species richness per site and jaccard indices between
# sites

## ---- load ----
PAstation <- readRDS("intermediate/02_species_diversity/PA_Mat_GaspObis_station.RDS")
PAsite <- readRDS("intermediate/02_species_diversity/PA_Mat_GaspObis_site.RDS")
PAall <- readRDS("intermediate/02_species_diversity/PA_Mat_GaspObis_allstations.RDS")
list_communities <- readRDS("intermediate/02_species_diversity/List_community.RDS")


## ---- initiate ----
list_sd_alpha_site <- list()
list_sd_alpha_station <- list()
list_sd_gamma <- list()
list_sd_beta_site_pair <- list()
list_sd_beta_site_multi <- list()
list_sd_beta_station_pair <- list()
list_sd_beta_station_multi <- list()


for (comm in names(list_communities)){
  
  print(comm)

  ## ---- subset to community ----
  PAstation_comm <-
    PAstation %>% 
    select(all_of(c("station", list_communities[[comm]])))
  
  PAsite_comm <-
    PAsite %>% 
    select(all_of(c("site", list_communities[[comm]])))
  
  PAall_comm <-
    PAall %>% 
    select(all_of(list_communities[[comm]]))
  
  ## ---- alpha diversity ----
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
  
  ## ---- gamma diversity ----
  list_sd_gamma[[comm]] <-
    data.frame(community = comm) %>% 
    cbind(richness_allstations = rowSums(PAall_comm))
    

  ## ---- beta diversity ----
  # beta diversity station level
  list_sd_beta_station_multi[[comm]] <- 
    beta.multi(PAstation_comm[,-1], index.family="jaccard")
  
  list_sd_beta_station_pair[[comm]] <- 
    beta.pair(PAstation_comm[,-1], index.family="jaccard")
  
  # beta diversity site level
  list_sd_beta_site_multi[[comm]] <- 
    beta.multi(PAsite_comm[,-1], index.family="jaccard")
  
  list_sd_beta_site_pair[[comm]] <- 
    beta.pair(PAsite_comm[,-1], index.family="jaccard")
  
}



## ---- collate ----
# alpha station
table_sd_alpha_station <- 
  do.call(rbind.data.frame, list_sd_alpha_station)

# alpha site
table_sd_alpha_site <- 
  do.call(rbind.data.frame, list_sd_alpha_site)

# alpha mean station
table_sd_alpha_mean_station <- 
  table_sd_alpha_station %>%
  group_by(community) %>% 
  summarise(richness_station_mean = mean(richness_station), 
            .groups = "keep")

# alpha mean site
table_sd_alpha_mean_site <- 
  table_sd_alpha_site %>%
  group_by(community) %>% 
  summarise(richness_site_mean = mean(richness_site), 
            .groups = "keep")

# gamma
table_sd_gamma <- 
  do.call(rbind.data.frame, list_sd_gamma) %>% 
  as_tibble()


# beta
table_sd_beta_station <- do.call(rbind.data.frame, list_sd_beta_station_multi) 
colnames(table_sd_beta_station) <- 
  paste0(colnames(table_sd_beta_station), "_station")
table_sd_beta_station <-
  table_sd_beta_station %>% 
  rownames_to_column(var = "community")

table_sd_beta_site <- do.call(rbind.data.frame, list_sd_beta_site_multi) 
colnames(table_sd_beta_site) <- 
  paste0(colnames(table_sd_beta_site), "_site")
table_sd_beta_site <-
  table_sd_beta_site %>% 
  rownames_to_column(var = "community")



## ---- merge ----
table_sd_global <-
  table_sd_gamma %>% 
  merge(table_sd_alpha_mean_site) %>% 
  merge(table_sd_alpha_mean_station) %>% 
  merge(table_sd_beta_site) %>% 
  merge(table_sd_beta_station)



## ---- export ----
write.csv(table_sd_global, "results/02_species_diversity/sd_table_global.csv", row.names = FALSE)
write.csv(table_sd_alpha_site, "results/02_species_diversity/sd_table_site.csv", row.names = FALSE)
write.csv(table_sd_alpha_station, "results/02_species_diversity/sd_table_station.csv", row.names = FALSE)
saveRDS(list_sd_beta_site_pair, "results/02_species_diversity/sd_list_pairwise_site.RDS")
saveRDS(list_sd_beta_station_pair, "results/02_species_diversity/sd_list_pairwise_station.RDS")





