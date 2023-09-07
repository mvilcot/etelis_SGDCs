
# ---- parameters ----
level = "site"

# communities delineation
# comm_delin = "taxonomic_scale_datasp2"
comm_delin = "taxonomic_scale_Fishbase"
# comm_delin = "depth_category"
# comm_delin = "phylogenetic_distance"
list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
comm = names(list_communities)[2]
metricSD = "richness_site"


# ---- load ----
level = "site"

sd_alpha <- read.csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))
sd_global <- read_csv(paste0("results/2_species_diversity/sd_table_global_", comm_delin, ".csv"))
sd_beta <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_site_", comm_delin, ".RDS"))

dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".RDS"))


# ---- analyse ----
temp <-
  dist_merge %>% 
  dplyr::select(contains(c("beta.jtu", "beta.jac", "beta.jne")))


temp <- 
  dist_merge %>% 
  dplyr::filter(!grepl("Hawaii", site))

summary(temp)

# order_sites <- 
#   data_sites[[level]] %>% 
#   levels(unique(data_sites[order(data_sites$order),][[level]])) %>% 
#   droplevels()
# data_samples$site




