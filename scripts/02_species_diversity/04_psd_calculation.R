# calculate phylogenetic metrics on inter-species phylogenies and output as a
# single table for oceans and sites

# load ----
phylo <-
  readRDS("./intermediate/sd_phylogenies.rds")

sd_table_global <-
  read_csv("./intermediate/sd_table_global.csv",
           show_col_types = FALSE)

sd_table_site <-
  read_csv("./intermediate/sd_table_site.csv",
           show_col_types = FALSE)

# calculate ----


# loop through locations
for(loc in names(phylo_fam)){
  
  # target phylogeny
  phylo_location <- phylo_fam[[loc]]
  
  # skip if there is only one species
  if(length(phylo_location$tip.label) > 1){
    # create distance matrix
    dist_location <-
      cophenetic.phylo(phylo_location)
    dist_location <-
      dist_location[upper.tri(dist_location)]
    
    # create psuedo-community matrix (all species present)
    com_mat <-
      matrix(rep(1,length(phylo_location$tip.label)),
             nrow = 1)
    colnames(com_mat) <- phylo_location$tip.label
    rownames(com_mat) <- c("global")
    
    # midpoint root the tree
    phylo_location <-
      midpoint.root(phylo_location)
    
    # pd
    pd <- 
      c(pd,
        as.numeric(pd(com_mat,phylo_location)[1]))
    
    # mpd
    mpd <-
      c(mpd,
        mean(dist_location))
    
    # vpd
    vpd <-
      c(vpd,
        var(dist_location))
    
  } else {
    # pd
    pd <- 
      c(pd,NA)
    
    # mpd
    mpd <-
      c(mpd,NA)
    
    # vpd
    vpd <-
      c(vpd,NA)
  }
  
  # meta
  meta_scale <-
    c(meta_scale,
      scale)
  meta_family <-
    c(meta_family,
      fam)
  meta_location <-
    c(meta_location,
      loc)
  
} # location loop


# compile ----
# put into single data frame
phylo_metrics <-
  data.frame(scale    = meta_scale,
             family   = meta_family,
             location = meta_location,
             sd_pd    = pd,
             sd_mpd   = mpd,
             sd_vpd   = vpd)
# ocean
sd_table_global$dummy <-
  paste0(sd_table_global$family,sd_table_global$ocean)

phylo_ocean <-
  phylo_metrics %>% 
  filter(scale == "ocean") %>% 
  mutate(dummy = paste0(family,location)) %>% 
  dplyr::select(-c(scale,family,location))

sd_table_global <-
  sd_table_global %>% 
  left_join(phylo_ocean, by = "dummy") %>% 
  dplyr::select(-dummy)

# site
sd_table_site$dummy <-
  paste0(sd_table_site$family,sd_table_site$sample_site)

phylo_site <-
  phylo_metrics %>% 
  filter(scale == "site") %>% 
  mutate(dummy = paste0(family,location)) %>% 
  dplyr::select(-c(scale,family,location))

sd_table_site <-
  sd_table_site %>% 
  left_join(phylo_site, by = "dummy") %>% 
  dplyr::select(-dummy)


# export ----
write_csv(sd_table_global,
          "./intermediate/sd_table_global.csv")

write_csv(sd_table_site,
          "./intermediate/sd_table_site.csv")