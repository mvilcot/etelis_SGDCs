library(tidyterra)


# ---- load ----
coord_vect <- vect(coord_site, geom=c("longitude", "latitude"), crs = "EPSG:4326")
MEOW <- vect("data/seascape/MEOW/Marine_Ecoregions_Of_the_World__MEOW_.shp", crs = "EPSG:3857")
# from https://geospatial.tnc.org/datasets/ed2be4cf8b7a451f84fd093c2e7660e3_0/explore
reefs <- vect("data/seascape/14_001_WCMC008_CoralReefs2021_v4_1/01_Data/WCMC008_CoralReef2021_Py_v4_1.shp")
# from https://data.unep-wcmc.org/datasets/1

# ---- transform ----
MEOW <- 
  MEOW %>% 
  project("EPSG:4326")

# ---- plot MEOW ----
ggplot() +
  geom_spatvector(data = MEOW) +
  geom_spatvector(data = coord_vect, color = "red") +
  theme_light() +
  crs

# ---- extract region ----
temp <- terra::extract(MEOW, coord_vect)
coord_site <- cbind(coord_site, temp)
coord_site


# ---- compute reef area by region ----
coord_site$area <- NA
gg_list <- list()
reef_list <- list()

for(i in 1:nrow(coord_site)){
  MEOWsub <- 
    terra::subset(MEOW, MEOW$ECO_CODE == coord_site$ECO_CODE[i])
  
  reef_list[[i]] <- terra::crop(reefs, MEOWsub)
  
  # plot(reef)
  # plot(coord_vect[i], add = T, col = "red", cex = 5)
  
  coord_site$area[i] <- sum(terra::expanse(reef_list[[i]]))
  
  gg_list[[i]] <- 
    ggplot() +
    geom_spatvector(data = MEOWsub) +
    geom_spatvector(data = reef_list[[i]]) +
    geom_spatvector(data = coord_vect[i,], 
                    color = "orange", size = 5, alpha = 0.8) +
    ggtitle(paste(rownames(coord_site)[i], coord_site$ECOREGION[i], coord_site$ECO_CODE[i], sep = ", ")) +
    theme_light()
}

gg_grob <- arrangeGrob(grobs = gg_list, ncol=5)
plot(gg_grob)
ggsave(gg_grob, width = 20, height = 10, 
       filename = paste0("results/3_distance_metrics/reef_shapes.png"))



# tidyterra::autoplot(reef)



# ---- area from Robuchon et al 2019 ----
area <- data.frame(area = env$Surf_area_km) 
rownames(area) <- env$BASIN
areadhm <- matrix(nrow = nrow(area), ncol = nrow(area), dimnames = list(rownames(area), rownames(area)))
for (i in rownames(area))
{
  for (j in rownames(area))
  {
    N1 <- area$area[rownames(area)==i]
    N2 <-  area$area[rownames(area)==j]
    areadhm[i, j] <- (2*N1*N2)/(N1+N2)
  }
} 
areadhm <- as.data.frame(areadhm) # dhm metric from Serrouya et al. 2012


## dhm metric
areadhm <- as.dist(areadhm)
areadhm <- (areadhm-mean(areadhm))/sd(areadhm)
areadhm <- data.frame(as.matrix(areadhm))
