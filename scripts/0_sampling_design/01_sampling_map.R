
# ---- setup sites data ----

## get bathymetry data
Bathy <- 
  getNOAA.bathy(lon1 = -180, lon2 = 180,
                lat1 = -90, lat2 = 90,
                resolution = 2,
                antimeridian = TRUE, # to be centered around pacific
                keep = TRUE)

## get number of Etelis coruscans sample by site
temp <- as.data.frame(table(data_samples$site))
colnames(temp) <- c("site", "N")

## merge to coordinates
coord_siteN <-
  coord_site %>% 
  rownames_to_column("site") %>% 
  left_join(temp) %>% 
  as_tibble()
write_csv(coord_siteN, "intermediate/0_sampling_design/coord_sites_N.csv")

## relevel sites
coord_siteN$site <- factor(coord_siteN$site, levels = coord_siteN$site)


# ---- data range Etelis coruscans ----

## Keep only presence data
data_EtelisP <- data_Etelis[data_Etelis$Etelis_coruscans == 1, ]



# ---- mapview ----

## transform to sf ----
coord_site_sf <- st_as_sf(coord_siteN, 
                        coords = c("longitude", "latitude"), 
                        crs = "WGS84", remove = FALSE)

data_sites_sf <- st_as_sf(data_sites, 
                          coords = c("Longitude_approx", "Latitude_approx"), 
                          crs = "WGS84", remove = FALSE)

## plot ----
mapsites <-
  mapview(shift.lonst(data_sites_sf), col.regions=list("grey"), cex = 2, legend = F) +
  mapview(shift.lonst(coord_site_sf), zcol = "site", cex = "N", label = "N")
mapsites

mapshot(mapsites, url = "results/0_sampling_design/Mapview_sampling_sites.html")


# ---- ggplot map ----
## load maps ----
# full map
map1 <-
  fortify(maps::map(fill=TRUE, plot=FALSE)) %>% 
  as_tibble()

# add map +360
map2 <- 
  map1 %>% 
  mutate(long = long + 360,
         group = group + max(map1$group) + 1)

## crop lon & lat extent ----
map <- 
  rbind(map1, map2) %>% 
  filter(long > 20  & long < 210 & lat <48 & lat > -48)



## plot ----
ggplot() +
  geom_point(data = shift.lon(data_EtelisP), aes(x = Longitude, y = Latitude),
             color = "grey", shape = 15, size = 1.8) +
  geom_polygon(data = map, aes(x = long, y = lat, group = group)) +
  geom_point(data = shift.lonst(coord_site_sf), 
             aes(x = longitude, y = latitude, color = site, size = N), 
             fill = "white", shape = 21) +
  ggrepel::geom_text_repel(data = shift.lonst(coord_site_sf),
            aes(x = longitude, y = latitude, color = site, label = site),
            hjust=0, vjust=0, max.overlaps = 10) +
  scale_color_viridis(discrete = T) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  coord_map() 

ggsave("results/0_sampling_design/ggplot_Map_sampling_sites.png", height = 8, width = 12)



# # ---- *** DRAFTS ----
# ## *** ggplot rnaturalearth ----
# world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
# worldshift <- st_shift_longitude(world)
# ggplot(data = worldshift) +
#   geom_sf() +
#   geom_point(data = shift.lonst(coord_site_sf), aes(x = longitude, y = latitude, color=site), size=1) +
#   scale_color_viridis(discrete = T) +
#   coord_sf(xlim = c(50, 250),
#            ylim = c(-75, 75)) 
# 
# 
# ## *** ggplot bathy rnarutalearth not shifted ----
# library(rnaturalearth)
# 
# # Get bathymetric data
# bat <- getNOAA.bathy(-12, -5, 35, 44, res = 4, keep = TRUE)
# bat_xyz <- as.xyz(bat)
# 
# # Import country data
# country <- ne_countries(scale = "medium", returnclass = "sf")
# 
# # Plot using ggplot and sf
# ggplot() + 
#   geom_sf(data = country) +
#   geom_tile(data = bat_xyz, aes(x = V1, y = V2, fill = V3)) +
#   # geom_contour(data = bat_xyz,
#   #              aes(x = V1, y = V2, z = V3),
#   #              breaks = -200, color = "grey85", linewidth = 0.5) +
#   geom_sf(data = country) +
#   geom_point(data = (coord_siteN), aes(x = longitude, y = latitude, colour = site)) +
#   scale_color_viridis(discrete = T) +
#   labs(x = "Longitude", y = "Latitude", fill = "Depth (m)") +
#   theme_minimal()
# 
# 
# ## *** auto bathy shifted ----
# autoplot.bathy(Bathy, geom=c("tile"), coast=T) +
#   # scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen") +
#   geom_point(data = shift.lon(coord_siteN), aes(x = longitude, y = latitude, colour = site)) +
#   labs(y = "Latitude", x = "Longitude", fill = "Elevation") 
# 

