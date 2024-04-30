
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


