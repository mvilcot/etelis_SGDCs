level = "site"
comm_delin = "taxonomic_scale"

gd_beta <- readRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))
sd_beta <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_", level, "_", comm_delin, ".RDS"))


## ---- GD ----
list_GDbeta <- list()

for (metricGD in names(gd_beta)){
  
  # get distance matrix
  mat_GDbeta <- as.matrix(gd_beta[[metricGD]])
  
  # order rows alphabetically
  mat_GDbeta <- mat_GDbeta[order(rownames(mat_GDbeta)), order(colnames(mat_GDbeta))]
  
  # pivot longer distance matrix
  melt_GDbeta <- melt.dist(dist = mat_GDbeta, metric = metricGD)
  
  # put into list
  list_GDbeta[[metricGD]] <- melt_GDbeta
  
}

# merge two distance matrix into one df
merge_gd_beta <- 
  plyr::join_all(list_GDbeta,
                 by = c(paste0(level, "1"), paste0(level, "2")))


## ---- SD ----
list_SDbeta <- list()

for (comm in names(sd_beta)){
  for (metricSD in names(sd_beta[[comm]])){
    
    # get distance matrix
    mat_SDbeta <- as.matrix(sd_beta[[comm]][[metricSD]])
    
    # order rows alphabetically
    mat_SDbeta <- mat_SDbeta[order(rownames(mat_SDbeta)), order(colnames(mat_SDbeta))]
    
    # pivot longer distance matrix
    melt_SDbeta <- melt.dist(dist = mat_SDbeta, metric = paste0(comm, '-', metricSD))
    
    # put into list
    list_SDbeta[[paste0(comm, '-', metricSD)]] <- melt_SDbeta
  }
}

# merge two distance matrix into one df
merge_sd_beta <- 
  plyr::join_all(list_SDbeta,
                 by = c(paste0(level, "1"), paste0(level, "2")))


## ---- merge ----

merge_beta <-
  full_join(merge_gd_beta,
            merge_sd_beta)

