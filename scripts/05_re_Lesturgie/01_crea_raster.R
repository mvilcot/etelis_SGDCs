#--------------------------------------------------#
#---------------- Creating raster -----------------#
#--------------------------------------------------#


## ---- COORDINATES + THETA LOCATIONS ---- 
# get mean coordinates by location
gd_alpha <- read.csv("results/01_genetic_diversity/gd_table_site.csv")

# get coord by site
coord_sites <- 
  data_sites %>% 
  # dplyr::select(all_of(level), Longitude_approx, Latitude_approx) %>% 
  # dplyr::rename("longitude" = "Longitude_approx",
  #        "latitude" = "Latitude_approx") %>%
  group_by(site) %>% 
  summarise(longitude = mean(Longitude_approx),
            latitude = mean(Latitude_approx)) %>% 
  left_join(gd_alpha) %>% 
  drop_na(Hs) %>% 
  column_to_rownames("site") %>% 
  as.data.frame()

# transform longitude to Pacific centered
lon <- coord_sites$longitude
coord_sites$longitude <- ifelse(lon < 0, lon + 360, lon)




## ---- empty raster ----
raster_empty <- raster(res = 0.6)


## ---- coord of reefs: convert in a raster ----
# reefs coordinates provided by ReefBase (http://www.reefbase.org)
coord_reefs <- read.csv("data/Lesturgie_2023/raster/Reefs/ReefLocations.csv")
coord_reefs <- cbind(coord_reefs$LON, coord_reefs$LAT)
raster_corals <- rasterize(coord_reefs, raster_empty, field=1)

# area_corals <- area(raster_corals)
# area_cells_c <- sqrt(mean(area_corals[]))

# attribute reefs value to 2000
A2 <- raster_corals@data@values
A2[which(is.na(A2)==F)] <- 2000
raster_corals <- setValues(raster_corals, A2)


## ---- coord of seamounts: convert in a raster ----
# downloaded from https://data.unep-wcmc.org/datasets/41
shp_seamounts <- shapefile("data/Lesturgie_2023/raster/Seamounts/Seamounts.shp") 
coord_seamounts <- cbind(shp_seamounts$LONG, shp_seamounts$LAT)
raster_seamounts <- rasterize(coord_seamounts, raster_empty, field=1)

# attribute seamounts value to 1000
A1 <- raster_seamounts@data@values
A1[which(is.na(A1)==F)] <- 1000
raster_seamounts <- setValues(raster_seamounts, A1)


## ---- loading a raster of the world ----
# downloaded at https://www.naturalearthdata.com/downloads/50m-raster-data/50m-natural-earth-1/
raster_world <- raster("data/Lesturgie_2023/raster/NE1_50M_SR_W/NE1_50M_SR_W.tif")

# area_world <- area(raster_world)
# area_cells_w <- sqrt(mean(area_world[]))
# fa <- area_world/area_cells_w

raster_world <- aggregate(raster_world, fact=18, fun=mean)


## ---- merge 3 rasters + pacific centering ----
raster_merge_0 <- merge(raster_corals, raster_seamounts, raster_world)
raster_merge_W <- crop(raster_merge_0, extent(-180, 0, -90, 90))
raster_merge_E <- crop(raster_merge_0, extent(0, 180, -90, 90))   
extent(raster_merge_W) <- c(180, 360, -90, 90)

raster_merge_PA <- merge(raster_merge_W, raster_merge_E)



## ---- just put different values to ocean and land ----
v <- raster_merge_PA$layer@data@values
w <- v
for (i in 1:length(v)){
  if (v[i]<1000){   # if neither coral or seamounts
    if (v[i]<=140){w[i] <- 500} # 500 for ocean (<140m)
    if (v[i]>140){w[i] <- 0.000001} # 0.000001 for land 
  }
}

raster_final <- setValues(raster_merge_PA, w)

# pdf("map_seamounts_corals_ocean.pdf")
# plot(raster_final, col=c("black","lightblue2","red4","red4","yellow"))
# points(coord_sites[,1:2], pch=16, cex=.5)
# legend(x = "bottom", col = c("black","lightblue2","red4","yellow"),
#       legend = c("land","ocean","seamounts","corals"), pch=16, ncol = 2)
# dev.off()



## ---- ggplot ----

df_final <- as.data.frame(raster_final, xy = TRUE)
cat <- df_final$layer
cat[cat==0.000001] <- "land"
cat[cat==500] <- "ocean"
cat[cat==1000] <- "seamount"
cat[cat==2000] <- "reef"

df_final$category <- cat


# ggplot() +
#   geom_raster(data = df_final, aes(x = x, y = y, fill = category)) +
#   scale_fill_manual(values = c("#233143","#E0EAEB","#EDAE49","#d1495b")) +
#   # "black","lightblue2","red4","yellow"
#   # scale_fill_viridis_d() +
#   geom_point(data = coord_sites, aes(longitude, latitude), color = "darkgreen") +
#   theme_void() +
#   coord_fixed()
# ggsave("results/05_re_Lesturgie/map_seamounts_corals_ocean.pdf",
#        height = 8, width = 16)


## ---- Let's crop to have only IndoP ----
raster_final_crop <- crop(raster_final, extent(35, 280, -35, 35))

cat("Construction raster complete","\n")

