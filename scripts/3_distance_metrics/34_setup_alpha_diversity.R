
# ---- parameters ----
level = "site"

# communtity delineation
comm_delin = "taxonomy"


list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
metricSD = "richness_site"
metricGD = "Hs"


# ---- load ----
gd_alpha <- read_csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"))
sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))

sd_alpha <-
  sd_alpha %>% 
  tidyr::pivot_wider(names_from = community, 
                     values_from = richness_site)

alpha_merge <-
  full_join(gd_alpha, sd_alpha, by = level)


## Rename Eupercaria ----
colnames(alpha_merge) <- gsub("Eupercaria/misc", "Eupercaria", colnames(alpha_merge))


# ---- save ----
alpha_merge %>% 
  write_csv(paste0("results/3_distance_metrics/alpha_diversity_", comm_delin, ".csv"))
