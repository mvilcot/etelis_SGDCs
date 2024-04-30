
# ---- load data ----
# parameters
level = "site"

# communtity delineation
comm_delin = "taxonomy"


# read
dist_mat <- readRDS("intermediate/3_distance_metrics/dist_geo_envtbdmean.RDS")
gd_mat <- readRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))
sd_mat <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_", level, "_", comm_delin, ".RDS"))

# unlist sd
sd_mat <- unlist(sd_mat, recursive = FALSE)


# ---- list of distance matrix ----
# add to one list
dist_mat <- 
  dist_mat %>% 
  append(gd_mat) %>% 
  append(sd_mat)


# ---- melt ----
dist_melt <- list()
for (i in 1:length(dist_mat)){
  dist_melt[[i]] <- melt.dist(distmat = dist_mat[[i]], metric = names(dist_mat)[[i]])
}

# ---- merge into a table ----
# full join
dist_merge <- 
  dist_melt %>%
  reduce(full_join, by = c("site1", "site2")) 

# create one column for both locations info
dist_merge$site <- 
  paste(dist_merge[["site1"]], 
        dist_merge[["site2"]],
        sep = "-")

# relocate to first column
dist_merge <- 
  dist_merge %>% 
  relocate("site")


## Rename Eupercaria!!!!!!!! ----
colnames(dist_merge) <- gsub("Eupercaria/misc", "Eupercaria", colnames(dist_merge))
names(dist_mat) <- gsub("Eupercaria/misc", "Eupercaria", names(dist_mat))



# ---- save ----
dist_merge %>% 
  write_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_mat %>% 
  saveRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))

