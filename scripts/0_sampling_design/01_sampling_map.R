
# data range Etelis coruscans ----
# Keep only presence data
data_EtelisP <- data_Etelis[data_Etelis$Etelis_coruscans == 1, ]


# transform to sf ----
vect_sites <- 
  data_sites %>% 
  vect(geom=c("longitude", "latitude"), crs = "WGS84") %>% 
  terra::rotate(left = FALSE)

vect_stations <- 
  data_stations %>% 
  vect(geom = c("Longitude_approx", "Latitude_approx"), crs = "WGS84") %>% 
  terra::rotate(left = FALSE)




# ---- Mapview sampling sites ----

## plot ----
mapsites <-
  mapview(vect_stations, col.regions=list("grey"), cex = 2, legend = F) +
  mapview(vect_sites, zcol = "site", cex = "number_samples", label = "N")
mapsites

mapshot(mapsites, url = "results/0_sampling_design/Mapview_sampling_sites.html")


# ---- Map bathymetry ----
#load
Bathy <- 
  getNOAA.bathy(lon1 = -180, lon2 = 180,
                lat1 = -50, lat2 = 50,
                resolution = 10,
                antimeridian = TRUE, # to be centered around pacific
                keep = TRUE)

Bathy %>% saveRDS("intermediate/0_sampling_design/bathymetry_10.RDS")
Bathy <- readRDS("intermediate/0_sampling_design/bathymetry_10.RDS")

# save bathymetry map
color_blues <- 
  colorRampPalette(c("purple", "blue", "cadetblue1", "cadetblue2", "white"))

pdf(file = "results/0_sampling_design/Bathymetry_map_10.pdf", height=5, width=12.7)
# png(file = "results/0_sampling_design/Bathymetry_map_10.png",
#     height=5, width=12.7, units="in", res=300)
plot.bathy(Bathy, 
           step=50000,
           deepest.isobath = -50000, shallowest.isobath = 0, 
           image = TRUE, bpal = color_blues(20))

scaleBathy(Bathy, 
           deg = 20, x = "bottomleft", inset = 2)

points(shift.lon(data_sites)$longitude, shift.lon(data_sites)$latitude,
       pch = 21, col = "black", bg = "red", cex = 1)

dev.off()




# ---- Map sampling sites ----

## load maps ----
# full map
map1 <-
  fortify(maps::map(fill=TRUE, plot=FALSE)) %>% 
  as_tibble()

# add map +360
## >>> SEE IF POSSIBLE TO DO WITH terra::rotate(left = FALSE) #########
map2 <- 
  map1 %>% 
  mutate(long = long + 360,
         group = group + max(map1$group) + 1)

# crop lon & lat extent
map <- 
  rbind(map1, map2) %>% 
  filter(long > 20  & long < 210 & lat <48 & lat > -48)

# tranform bathymetry data
bat_xyz <- as.xyz(Bathy)


## plot ----
# library(sf)
# autoplot.bathy(Bathy, geom="tile")

## bathy
gg <- 
  ggplot() +
  
  ## Bathymetry
  geom_tile(data = bat_xyz, aes(x = V1, y = V2, fill = V3),
            width = 2) +
  scale_fill_gradient2(low="dodgerblue4", high="white") +
  
  ## Etelis presence
  geom_point(data = shift.lon(data_EtelisP), aes(x = Longitude, y = Latitude),
             color = "grey40", shape = 15, size = 1.8) +
  
  ## Countries
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
  
  ## Sites
  geom_spatvector(data = vect_sites, aes(color = site, size = number_samples),
                  fill = "white", shape = 21, show.legend = T) +
  # geom_spatvector_text(data = vect_sites, aes(label = site, color = site)) +
  ggrepel::geom_text_repel(data = shift.lon(data_sites),
            aes(x = longitude, y = latitude, color = site, label = site),
            hjust=0, vjust=0, max.overlaps = 10) +
  scale_color_viridis(discrete = T) +
  
  ## Theme
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(30, 200), ylim = c(-45, 40),
           crs = "EPSG:3832")

gg

ggsave("results/0_sampling_design/ggplot_Map_sampling_sites_bathy_10.png", 
       gg,
       height = 8, width = 12)


## no bathy
# gg <- 
#   ggplot() +
#  
#   ## Etelis presence
#   geom_point(data = shift.lon(data_EtelisP), aes(x = Longitude, y = Latitude),
#              color = "grey40", shape = 15, size = 1.8) +
#   
#   ## Countries
#   geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
#   
#   ## Sites
#   geom_spatvector(data = vect_sites, aes(color = site, size = number_samples),
#                   fill = "white", shape = 21, show.legend = T) +
#   # geom_spatvector_text(data = vect_sites, aes(label = site, color = site)) +
#   ggrepel::geom_text_repel(data = shift.lon(data_sites),
#                            aes(x = longitude, y = latitude, color = site, label = site),
#                            hjust=0, vjust=0, max.overlaps = 10) +
#   scale_color_viridis(discrete = T) +
#   
#   ## Theme
#   theme_minimal() +
#   theme(legend.position = "none") +
#   labs(x = "Longitude", y = "Latitude")
# gg

gg <- 
  ggplot() +
  
  ## Etelis presence
  geom_point(data = shift.lon(data_EtelisP), aes(x = Longitude, y = Latitude),
             color = "darkgrey", shape = 15, size = 1.8) +
  ## Countries
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
  ## Sites
  geom_point(data = shift.lon(data_sites),
             aes(x = longitude, y = latitude, color = site, size = number_samples),
             fill = "white", shape = 21) +
  ggrepel::geom_text_repel(data = shift.lon(data_sites),
                           aes(x = longitude, y = latitude, color = site, label = site),
                           hjust=0, vjust=0, max.overlaps = 5, nudge_x = -10) +
  ## Theme
  scale_color_viridis(discrete = T) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  coord_map

gg

ggsave("results/0_sampling_design/ggplot_Map_sampling_sites.png", 
       gg,
       height = 7, width = 11)


## no range
#### >>> HERE ####
reefs <- vect("data/seascape/14_001_WCMC008_CoralReefs2021_v4_1/01_Data/WCMC008_CoralReef2021_Py_v4_1.shp")
# from https://data.unep-wcmc.org/datasets/1
reefsC <-
  reefs %>% 
  terra::rotate(left = FALSE) %>% 
  terra::crop(extent(20, 210, -48, 48))


gg <- 
  ggplot() +
  
  ## Countries
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
   
  ## reefs
  geom_spatvector(data = reefsC, color = "grey40", linewidth = 0.02) +
  
  ## Sites
  geom_point(data = shift.lon(data_sites),
             aes(x = longitude, y = latitude, fill = site, color = site,size = number_samples),
             alpha = 0.6, shape = 21) +
  ggrepel::geom_text_repel(data = shift.lon(data_sites),
                           size = 3.2,
                           aes(x = longitude, y = latitude, color = site, label = site),
                           hjust=0.5, vjust=0, max.overlaps = 10, nudge_x = -5, nudge_y = 0,
                           direction = "both",
                           force = 5,
                           bg.color = "grey70", bg.r = 0.02) +

  ## Theme
  scale_fill_viridis(discrete = T, end = 1) +
  scale_color_viridis(discrete = T, end = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(plot.background = element_rect(fill = 'white', color = "white")) +
  labs(x = "Longitude", y = "Latitude")
  # coord_cartesian()

gg

ggsave("results/0_sampling_design/ggplot_Map_sampling_sites_norange_REEFS_final.png", 
       gg, height = 7, width = 11, dpi = 500)
ggsave("results/0_sampling_design/ggplot_Map_sampling_sites_norange_REEFS_final.pdf", 
       gg, height = 7, width = 11)


# ---- Map seamounts, reefs ----
## empty raster ----
raster_empty <- raster::raster(res = 0.6)


## coord of reefs: convert in a raster ----
# reefs coordinates provided by ReefBase (http://www.reefbase.org)
coord_reefs <- read_csv("data/Lesturgie_2023/raster/Reefs/ReefLocations.csv")
coord_reefs <- cbind(coord_reefs$LON, coord_reefs$LAT)
raster_corals <- raster::rasterize(coord_reefs, raster_empty, field=1)

# attribute reefs value to 2000
A2 <- raster_corals@data@values
A2[which(is.na(A2)==F)] <- 2000
raster_corals <- raster::setValues(raster_corals, A2)


## coord of seamounts: convert in a raster ----
# downloaded from https://data.unep-wcmc.org/datasets/41
shp_seamounts <- shapefile("data/Lesturgie_2023/raster/Seamounts/Seamounts.shp") 
coord_seamounts <- cbind(shp_seamounts$LONG, shp_seamounts$LAT)
raster_seamounts <- raster::rasterize(coord_seamounts, raster_empty, field=1)

# attribute seamounts value to 1000
A1 <- raster_seamounts@data@values
A1[which(is.na(A1)==F)] <- 1000
raster_seamounts <- raster::setValues(raster_seamounts, A1)


## raster of the world ----
# downloaded at https://www.naturalearthdata.com/downloads/50m-raster-data/50m-natural-earth-1/
raster_world <- raster("data/Lesturgie_2023/raster/NE1_50M_SR_W/NE1_50M_SR_W.tif")
raster_world <- raster::aggregate(raster_world, fact=18, fun=mean)


## merge 3 rasters + pacific centering ----
raster_merge_0 <- raster::merge(raster_corals, raster_seamounts, raster_world)
raster_merge_W <- raster::crop(raster_merge_0, extent(-180, 0, -90, 90))
raster_merge_E <- raster::crop(raster_merge_0, extent(0, 180, -90, 90))   
extent(raster_merge_W) <- c(180, 360, -90, 90)

raster_merge_PA <- raster::merge(raster_merge_W, raster_merge_E)



## attribute values to ocean and land ----
raster1 <- raster_merge_PA$layer@data@values
raster2 <- raster1
for (i in 1:length(raster1)){
  if (raster1[i]<1000){   # if neither coral or seamounts
    if (raster1[i]<=140){raster2[i] <- 1} # 1 for ocean (<140m)
    if (raster1[i]>140){raster2[i] <- 0} # 0 for land 
  }
}

raster_final <- raster::setValues(raster_merge_PA, raster2)




## plot reefs + seamounts ----
df_final <- as.data.frame(raster_final, xy = TRUE)
df_final$category <- df_final$layer
df_final$category[df_final$category==0] <- "land"
df_final$category[df_final$category==1] <- "ocean"
df_final$category[df_final$category==1000] <- "seamount"
df_final$category[df_final$category==2000] <- "reef"

ggplot() +
  geom_raster(data = df_final, aes(x = x, y = y, fill = category)) +
  scale_fill_manual(values = c("white","#dae7ed","#EDAE49","#5e7d4a")) +
  geom_point(data = shift.lon(data_sites), aes(longitude, latitude), size = 3, color = "red") +
  theme_void() +
  coord_sf(xlim = c(30, 220),
           ylim = c(-50, 50))
ggsave("results/0_sampling_design/Map_seamounts_corals_ocean.png",
       height = 8, width = 16, units = "in", dpi = 300)


# ---- Map reef shapes ----
## load ----
MEOW <- vect("data/seascape/MEOW/Marine_Ecoregions_Of_the_World__MEOW_.shp", crs = "EPSG:3857")
# from https://geospatial.tnc.org/datasets/ed2be4cf8b7a451f84fd093c2e7660e3_0/explore
reefs <- vect("data/seascape/14_001_WCMC008_CoralReefs2021_v4_1/01_Data/WCMC008_CoralReef2021_Py_v4_1.shp")
# from https://data.unep-wcmc.org/datasets/1


## extract region of each sampling site ----
MEOW <- 
  MEOW %>% 
  terra::project("WGS84")

vect_sitesWGS <- 
  data_sites %>% 
  vect(geom=c("longitude", "latitude"), crs = "WGS84") # back to -180/180

data_regions <-
  data_sites %>% 
  cbind(terra::extract(MEOW, vect_sitesWGS))


## plot MEOW ----
ggplot() +
  geom_spatvector(data = MEOW) +
  geom_spatvector(data = vect_sitesWGS, color = "#EDAE19") +
  theme_light()


## project to Pacific ----
MEOWP <-
  MEOW %>%
  terra::project("EPSG:3832") 

reefsP <-
  reefs %>%
  terra::project("EPSG:3832")

vect_sitesP <-
  vect_sites %>%
  terra::project("EPSG:3832")

vect_stationsP <-
  vect_stations %>% 
  terra::project("EPSG:3832")




## plot by sampling site ----

gg_list <- list()
reef_list <- list()
# data_regions$area <- NA

for(i in 1:nrow(data_regions)){
  focalsite <- as.character(data_regions$site[i])
  
  MEOWsub <- 
    terra::subset(MEOWP, MEOWP$ECO_CODE == data_regions$ECO_CODE[i])
  
  reef_list[[focalsite]] <- terra::crop(reefsP, MEOWsub)
  
  stations_sub <- 
    vect_stationsP %>% 
    dplyr::filter(site == focalsite)
  
  # ## basic plot
  # plot(reef_list[[i]])
  # plot(vect_sitesP[i], add = T, col = "red", cex = 5)
  
  ## compute area
  # data_regions$area[i] <- sum(terra::expanse(reef_list[[i]]))
  
  gg_list[[focalsite]] <-
    ggplot() +
    geom_spatvector(data = MEOWsub, fill = "grey95", color = "grey40") +
    geom_spatvector(data = reef_list[[i]], color = "grey10") +
    geom_spatvector(data = stations_sub, shape = 3,
                    color = "cyan4", size = 3, alpha = 0.8) +
    geom_spatvector(data = vect_sitesP[i,], shape = 19,
                    fill = "cyan4", color = "cyan4", size = 4, alpha = 0.9) + #EDAE19
    # ggrepel::geom_text_repel(data = shift.lon(data_sites)[i,], aes(x=longitude, y=latitude, label = focalsite), color = "orange") +
    ggtitle(paste(focalsite, data_regions$ECOREGION[i], sep = ", ")) + 
    theme_light()
}


gg_grob <- arrangeGrob(grobs = gg_list, ncol=5)
plot(gg_grob)
ggsave(gg_grob, width = 20, height = 10, 
       filename = paste0("results/0_sampling_design/MEOW_reef_shapes_stations.png"))



## ---- *** DRAFTS ----
# 
# sf_sites <- st_as_sf(data_sites, 
#                           coords = c("longitude", "latitude"), 
#                           crs = "WGS84", remove = FALSE)
# 
# sf_stations <- st_as_sf(data_stations, 
#                           coords = c("Longitude_approx", "Latitude_approx"), 
#                           crs = "WGS84", remove = FALSE)
# 
### *** ggplot rnaturalearth ----
# world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
# worldshift <- st_shift_longitude(world)
# ggplot(data = worldshift) +
#   geom_sf() +
#   geom_point(data = shift.lon(sf_sites), aes(x = longitude, y = latitude, color=site), size=1) +
#   scale_color_viridis(discrete = T) +
#   coord_sf(xlim = c(50, 250),
#            ylim = c(-75, 75)) 
# 
# 
### *** ggplot bathy rnarutalearth not shifted ----
# library(rnaturalearth)
# 
# # Get bathymetric data
# Bathy <- readRDS("intermediate/0_sampling_design/bathymetry.RDS")
# bat_xyz <- as.xyz(Bathy)
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
#   geom_point(data = (data_sites), aes(x = longitude, y = latitude, colour = site)) +
#   scale_color_viridis(discrete = T) +
#   labs(x = "Longitude", y = "Latitude", fill = "Depth (m)") +
#   theme_minimal()
# 
# 
### *** auto bathy shifted ----
# autoplot.bathy(Bathy, geom=c("tile"), coast=T) +
#   # scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen") +
#   geom_point(data = shift.lon(data_sites), aes(x = longitude, y = latitude, colour = site)) +
#   labs(y = "Latitude", x = "Longitude", fill = "Elevation") 
# 

