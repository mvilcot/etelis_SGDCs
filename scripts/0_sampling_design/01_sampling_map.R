library(mapview)
library(ggspatial)

# ---- Prepare data ----
## reefs ----
# First download this file from the WCMC ressource https://data.unep-wcmc.org/datasets/1
reefs <- vect("data/seascape/14_001_WCMC008_CoralReefs2021_v4_1/01_Data/WCMC008_CoralReef2021_Py_v4_1.shp")

reefsC <-
  reefs %>% 
  terra::rotate(left = FALSE) %>% 
  terra::crop(extent(20, 210, -48, 48))


## maps ----
# full map
map1 <-
  fortify(maps::map(fill=TRUE, plot=FALSE)) %>% 
  as_tibble()

# add map +360
map2 <- 
  map1 %>% 
  mutate(long = long + 360,
         group = group + max(map1$group) + 1)

# crop lon & lat extent
map <- 
  rbind(map1, map2) %>% 
  filter(long > 20  & long < 210 & lat <48 & lat > -48)



# ---- Plot sampling sites ----
data_sites$labels <- paste0(data_sites$site, ' (', data_sites$number_samples, ')')
data_sites$labels <- gsub('_', ' ', data_sites$labels)
data_sites$labels <- gsub('W Australia', 'Western Australia', data_sites$labels)

gg <- 
  ggplot() +
  
  ## Countries
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
   
  ## reefs (optional if reef shp file not available)
  geom_spatvector(data = reefsC, color = "grey40", linewidth = 0.02) +
  
  ## Sites
  geom_point(data = shift.lon(data_sites),
             aes(x = longitude, y = latitude, fill = site, color = site), #size = number_samples
             size = 3, alpha = 0.6, shape = 21) +
  ggrepel::geom_text_repel(data = shift.lon(data_sites),
                           size = 3.1, segment.color = 'transparent',
                           aes(x = longitude, y = latitude, color = site, label = labels),
                           hjust = 0.9, vjust = 0, max.overlaps = 10, nudge_y = -3, 
                           direction = "both",
                           force = 1,
                           bg.color = "grey70", bg.r = 0.02) +

  ## Theme
  scale_fill_viridis(discrete = T, end = 1) +
  scale_color_viridis(discrete = T, end = 1) +
  theme_minimal() +
  theme(legend.position = 'none', 
        plot.background = element_rect(fill = 'white', color = "white")) +
  guides(color = 'none', fill = 'none') +
  labs(x = "", y = "") +

  # scales
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.4, "in"), width = unit(0.4, "in"),
                         pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_orienteering) 


#### {FIGURE 1} ####
ggsave("results/0_sampling_design/_1_ggplot_Map_sampling_sites_norange_REEFS_FINAL.png", 
       gg, height = 6, width = 11, dpi = 500)
ggsave("results/0_sampling_design/_1_ggplot_Map_sampling_sites_norange_REEFS_FINAL.pdf",
       gg, height = 6, width = 11)

