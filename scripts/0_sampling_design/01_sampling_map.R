# ---- Load ----
## setup sites data ----
# get number of Etelis coruscans sample by site
temp <- as.data.frame(table(data_samples$site))
colnames(temp) <- c("site", "N")

# merge to coordinates
coord_siteN <-
  coord_site %>% 
  rownames_to_column("site") %>% 
  left_join(temp) %>% 
  as_tibble()
write_csv(coord_siteN, "intermediate/0_sampling_design/coord_sites_N.csv")

# relevel sites
coord_siteN$site <- factor(coord_siteN$site, levels = coord_siteN$site)


##  data range Etelis coruscans ----
# Keep only presence data
data_EtelisP <- data_Etelis[data_Etelis$Etelis_coruscans == 1, ]



# ---- Mapview sampling sites ----

## transform to sf ----
coord_site_sf <- st_as_sf(coord_siteN, 
                        coords = c("longitude", "latitude"), 
                        crs = "WGS84", remove = FALSE)

data_sites_sf <- st_as_sf(data_sites, 
                          coords = c("Longitude_approx", "Latitude_approx"), 
                          crs = "WGS84", remove = FALSE)

## plot ----
mapsites <-
  mapview(shift.lon(data_sites_sf), col.regions=list("grey"), cex = 2, legend = F) +
  mapview(shift.lon(coord_site_sf), zcol = "site", cex = "N", label = "N")
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

points(shift.lon(coord_site)$longitude, shift.lon(coord_site)$latitude,
       pch = 21, col = "black", bg = "red", cex = 1)

dev.off()




# ---- Map sampling sites ----

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

# crop lon & lat extent
map <- 
  rbind(map1, map2) %>% 
  filter(long > 20  & long < 210 & lat <48 & lat > -48)

# tranform bathymetry data
bat_xyz <- as.xyz(Bathy)


## plot ----
library(sf)

## bathy
gg <- 
  ggplot() +
  # autoplot.bathy(Bathy, geom="tile")
  geom_tile(data = bat_xyz, aes(x = V1, y = V2, fill = V3),
            width = 2) +
  scale_fill_gradient2(low="dodgerblue4", high="white") +
  geom_point(data = shift.lon(data_EtelisP), aes(x = Longitude, y = Latitude),
             color = "darkgreen", shape = 15, size = 1.8) +
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
  geom_point(data = shift.lon(coord_site_sf),
             aes(x = longitude, y = latitude, color = site, size = N),
             fill = "white", shape = 21) +
  scale_color_viridis(discrete = T) +
  ggrepel::geom_text_repel(data = shift.lon(coord_site_sf),
            aes(x = longitude, y = latitude, color = site, label = site),
            hjust=0, vjust=0, max.overlaps = 10) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  coord_map() +
  coord_sf(xlim = c(30, 200), ylim = c(-45, 40))


ggsave("results/0_sampling_design/ggplot_Map_sampling_sites_bathy_10.png", 
       gg,
       height = 8, width = 12)


## no bathy
gg <- 
  ggplot() +
  geom_point(data = shift.lon(data_EtelisP), aes(x = Longitude, y = Latitude),
             color = "darkgrey", shape = 15, size = 1.8) +
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
  geom_point(data = shift.lon(coord_site_sf),
             aes(x = longitude, y = latitude, color = site, size = N),
             fill = "white", shape = 21) +
  ggrepel::geom_text_repel(data = shift.lon(coord_site_sf),
                           aes(x = longitude, y = latitude, color = site, label = site),
                           hjust=0, vjust=0, max.overlaps = 10) +
  scale_color_viridis(discrete = T) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  coord_map() 
  # coord_sf(xlim = c(30, 200), ylim = c(-45, 40))


ggsave("results/0_sampling_design/ggplot_Map_sampling_sites.png", 
       gg,
       height = 8, width = 12)




# ---- Map seamounts, reefs ----
## empty raster ----
raster_empty <- raster::raster(res = 0.2)


## coord of reefs: convert in a raster ----
# reefs coordinates provided by ReefBase (http://www.reefbase.org)
coord_reefs <- read.csv("data/Lesturgie_2023/raster/Reefs/ReefLocations.csv")
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
raster_world <- rnaturalearth::
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
  scale_fill_manual(values = c("white","#dae7ed","#EDAE49","#233143")) +
  geom_point(data = shift.lon(coord_site), aes(longitude, latitude), color = "red") +
  theme_void() +
  coord_fixed()
ggsave("results/0_sampling_design/Map_seamounts_corals_ocean.pdf",
       height = 8, width = 16)


## crop to have only IndoP ----
raster_final_crop <- raster::crop(raster_final, extent(35, 280, -35, 35))

cat("Construction raster complete \n")



# ---- Map reef shapes ----
library(tidyterra)

## load ----
coord_vect <- vect(coord_site, geom=c("longitude", "latitude"), crs = "EPSG:4326")
MEOW <- vect("data/seascape/MEOW/Marine_Ecoregions_Of_the_World__MEOW_.shp", crs = "EPSG:3857")
# from https://geospatial.tnc.org/datasets/ed2be4cf8b7a451f84fd093c2e7660e3_0/explore
reefs <- vect("data/seascape/14_001_WCMC008_CoralReefs2021_v4_1/01_Data/WCMC008_CoralReef2021_Py_v4_1.shp")
# from https://data.unep-wcmc.org/datasets/1

## project ----
MEOW <- 
  MEOW %>% 
  terra::project("EPSG:4326")

## extract region of each sampling site ----
temp <- terra::extract(MEOW, coord_vect)
coord_siteT <- cbind(coord_site, temp)
coord_siteT

## project to Pacific ----
MEOW <-
  MEOW %>%
  terra::project("EPSG:3832")

reefs <-
  reefs %>%
  terra::project("EPSG:3832")

coord_vect <-
  coord_vect %>%
  terra::project("EPSG:3832")

stations_vect <-
  vect(data_sites, geom=c("Longitude_approx", "Latitude_approx"), crs = "EPSG:4326") %>% 
  terra::project("EPSG:3832")




## plot MEOW ----
ggplot() +
  geom_spatvector(data = MEOW) +
  geom_spatvector(data = coord_vect, color = "orange") +
  theme_light() 



## plot by sampling site ----

gg_list <- list()
reef_list <- list()
# coord_siteT$area <- NA

for(i in 1:nrow(coord_siteT)){
  focalsite <- rownames(coord_siteT)[i]
  
  MEOWsub <- 
    terra::subset(MEOW, MEOW$ECO_CODE == coord_siteT$ECO_CODE[i])
  
  reef_list[[focalsite]] <- terra::crop(reefs, MEOWsub)
  
  stations_sub <- 
    stations_vect %>% 
    filter(site == focalsite)
  
  ## basic plot
  # plot(reef_list[[i]])
  # plot(coord_vect[i], add = T, col = "red", cex = 5)
  
  ## compute area
  # coord_siteT$area[i] <- sum(terra::expanse(reef_list[[i]]))
  
  gg_list[[focalsite]] <-
    ggplot() +
    geom_spatvector(data = MEOWsub, fill = "grey95", color = "darkgrey") +
    geom_spatvector(data = reef_list[[i]], color = "#EDAE49") +
    geom_spatvector(data = stations_sub,
                    color = "grey40", size = 2, alpha = 0.8) +
    geom_spatvector(data = coord_vect[i,],
                    color = "black", size = 5, alpha = 0.8) +
    # ggrepel::geom_text_repel(data = shift.lon(coord_site)[i,], aes(x=longitude, y=latitude, label = focalsite), color = "orange") +
    ggtitle(paste(rownames(coord_siteT)[i], coord_siteT$ECOREGION[i], sep = ", ")) + 
    theme_bw()
}

gg_grob <- arrangeGrob(grobs = gg_list, ncol=5)
# plot(gg_grob)
ggsave(gg_grob, width = 25, height = 12.5, 
       filename = paste0("results/0_sampling_design/MEOW_reef_shapes_shift_stations.png"))



# # ---- *** DRAFTS ----
# ## *** ggplot rnaturalearth ----
# world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
# worldshift <- st_shift_longitude(world)
# ggplot(data = worldshift) +
#   geom_sf() +
#   geom_point(data = shift.lon(coord_site_sf), aes(x = longitude, y = latitude, color=site), size=1) +
#   scale_color_viridis(discrete = T) +
#   coord_sf(xlim = c(50, 250),
#            ylim = c(-75, 75)) 
# 
# 
# ## *** ggplot bathy rnarutalearth not shifted ----
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

