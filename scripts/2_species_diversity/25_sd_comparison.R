
# ---- parameters ----
## communities delineation
comm_delin = "taxonomy"
# comm_delin = "taxonomy_depth1_crosses45-400m"
# comm_delin = "taxonomy_depth2_contains45-400m"
# comm_delin = "taxonomy_depth3_within45-400m"
# comm_delin = "taxonomy_env_reef-associated"
list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))
names_communities <- c("Etelinae", "Lutjanidae", "Eupercaria", "Teleostei")

# parameters
comm = names(list_communities)[2]
metricSD = "richness_site"


# ---- load ----
level = "site"

sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))
sd_global <- read_csv(paste0("results/2_species_diversity/sd_table_global_", comm_delin, ".csv"))
sd_beta <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_site_", comm_delin, ".RDS"))

dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".csv"))

dist_mat <-
  readRDS(paste0("intermediate/3_distance_metrics/dist_geo_envtbdmean_gd_sd_", comm_delin, ".RDS"))


# ---- summary ----
temp <-
  dist_merge %>% 
  dplyr::select(contains(c("beta.jtu", "beta.jac", "beta.jne")))

# temp <-
#   dist_merge %>%
#   dplyr::filter(!grepl("Hawaii", site))

summary(temp)

# order_sites <- 
#   data_stations[[level]] %>% 
#   levels(unique(data_stations[order(data_stations$order),][[level]])) %>% 
#   droplevels()
# data_samples$site


# ---- alpha boxplot ----
## richness station ----
# load
level = "station"
sd_alpha <- read_csv(paste0("results/2_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))

df <-
  sd_alpha %>% 
  left_join(data_stations, by = "station")

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
  data_sites %>% 
  left_join(sd_alpha, by = "site")

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
sd_alpha$community <- gsub('Eupercaria/misc', 'Eupercaria', sd_alpha$community)
df <-
  data_sites %>% 
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
  scale_fill_manual('', values = color_perso, labels = LABELS) +
  ylab('Species richness') +
  xlab('Longitude') +
  theme_light()

### {FIGURE S4} ####
ggsave(paste0("results/2_species_diversity/_S4_alpha_sd_richness_longitude_site.png"),
       width = 10, height = 3.5, dpi = 500)


## richness ~ dist to CT ----
### !!!!!!!! TO DO !!!!!!!!!!!!!!!! ----


# Violin plots by community ----
## alpha ----
sd_alpha <- 
  sd_alpha %>% 
  mutate(site = factor(site, levels = levels(data_sites$site))) %>% 
  mutate(community = gsub("/misc", "", community)) %>% 
  mutate(community = factor(community, levels = names_communities)) %>% 
  arrange(community, site)

sd_alpha %>% 
  ggplot(aes(community, richness_site)) +
  geom_violin(fill = "grey50", color = "grey50") +
  geom_boxplot(width=0.1) +
  # geom_hline(yintercept = 45, color = "orange") +
  xlab("") + ylab('α-diversity') +
  theme_light()

## beta ----
temp_beta <- 
  dist_merge %>% 
  dplyr::select(-c(Fst, GstPP.hed, D.Jost, jtu, jac, jne)) %>% 
  pivot_longer(cols = -c(site, site1, site2, geodist, seadist, environment), 
               names_to = "variable") %>% 
  separate(variable, c("community", "drop", "variable"), sep = "\\.") %>% 
  mutate(variable = paste0("β-", variable)) %>% 
  dplyr::filter(variable != "β-jne")

temp_beta %>% 
  ggplot(aes(community, value, fill = variable)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  # geom_hline(yintercept = 45, color = "orange") +
  xlab("") + ylab('β-diversity') +
  theme_light()


