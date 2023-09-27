
# ---- parameters ----
## communities delineation
comm_delin = "taxonomy"
# comm_delin = "taxonomy_depth1_crosses45-400m"
# comm_delin = "taxonomy_depth2_contains45-400m"
# comm_delin = "taxonomy_depth3_within45-400m"
# comm_delin = "taxonomy_env_reef-associated"
list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
comm = names(list_communities)[2]
metricSD = "richness_site"


# ---- load ----
level = "site"

sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))
sd_global <- read_csv(paste0("results/2_species_diversity/sd_table_global_", comm_delin, ".csv"))
sd_beta <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_site_", comm_delin, ".RDS"))

dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".RDS"))


# ---- summary ----
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


# ---- alpha boxplot ----
## richness station ----
# load
level = "station"
sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))

df <-
  sd_alpha %>% 
  left_join(data_sites)

df$community <- factor(df$community, levels = unique(df$community))

# plot
ggplot(df, aes(x=site, y=richness_station, color=site)) + 
  geom_boxplot() +
  facet_grid(rows = vars(community), scales ="free_y") +
  scale_color_manual(values = color_perso) +
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

# save
ggsave(paste0("results/2_species_diversity/boxplot_alpha_richness_station.png"),
       width = 8, height = 8)



## richness site ----
# load
level = "site"
sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))

df <-
  coord_site %>% 
  rownames_to_column("site") %>% 
  left_join(sd_alpha)

df$site <- factor(df$site, levels = unique(df$site))
df$community <- factor(df$community, levels = unique(df$community))
df <- shift.lon(df)

# plot
ggplot(df, aes(x=site, y=richness_site, color=site)) + 
  geom_boxplot() +
  facet_grid(rows = vars(community), scales ="free_y") +
  scale_color_manual(values = color_perso) +
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

# save
ggsave(paste0("results/2_species_diversity/boxplot_alpha_richness_site.png"),
       width = 8, height = 8)



## richness ~ longitude ----
# load
level = "site"
sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))

df <-
  coord_site %>% 
  rownames_to_column("site") %>% 
  left_join(sd_alpha)

df$site <- factor(df$site, levels = unique(df$site))
df$community <- factor(df$community, levels = unique(df$community))
df <- shift.lon(df)

# plot
ggplot(df, aes(x=longitude, y=richness_site, fill=site, group=site)) + 
  geom_bar(position = "dodge", stat="identity") +
  # geom_smooth(color="grey") +
  # geom_point(aes(color=site)) +
  facet_wrap(vars(community), ncol = 4, scales ="free_y") +
  scale_fill_manual(values = color_perso) +
  theme_light()

# save
ggsave(paste0("results/2_species_diversity/alpha_sd_richness_longitude_site.png"),
       width = 15, height = 6)


## richness ~ dist to CT ----
### !!!!!!!! TO DO !!!!!!!!!!!!!!!! ----



