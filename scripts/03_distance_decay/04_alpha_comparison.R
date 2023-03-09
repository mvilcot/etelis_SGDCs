

## ---- load data ----
level = "site"
comm_delin = "taxonomic_scale"

gd_alpha <- read.csv(paste0("results/01_genetic_diversity/gd_table_", level, ".csv"))
sd_alpha <- read.csv(paste0("results/02_species_diversity/sd_table_", level, "_", comm_delin, ".csv"))


## ---- handle data ----
# remove sites
gd_alpha <- 
  gd_alpha %>% 
  filter(site != "Seychelles")
  
sd_alpha <- 
  sd_alpha %>% 
  filter(site != "Seychelles")

# scale metrics
gd_alpha <- 
  gd_alpha %>%
  mutate_at(c("Ho", "Hs"), scale)

sd_alpha <- 
  sd_alpha %>%
  group_by(community) %>%
  mutate_at("richness_site", scale)


# merge tables
table_alpha <- 
  gd_alpha %>% 
  left_join(sd_alpha, by = level) %>% 
  select(site, community, Hs, richness_site) %>% 
  pivot_longer(cols = c("Hs", "richness_site"), names_to = "metric")


## ---- plot ----
ggplot(table_alpha, aes(x=site, y=value, color=metric)) + 
  geom_boxplot() +
  facet_wrap(~community, scale="free")







