# -----------  Comparison send samples to sequenced ones ------------------------

samples_sent <- read.csv("../Labwork/DART Shipping/dart_sample_trackingfile_template.csv")
samples_sequenced <- read.csv("Data/DartSeq/targets_H5MF3DMXY_2.csv")

samples_sent$sequenced <-  ifelse(samples_sent$Genotype %in% samples_sequenced$genotype, "yes", "NO")
# samples_sent$quality <- 
write.csv(samples_sent, "Samples_sequenced.csv", row.names = F)

  

# -----------  Map presence data ------------------------------------------------

library("ggplot2")
library("rnaturalearth")
library("rnaturalearthdata")
library("mapdata")
library("sf")
library("mapview")

## Read presence data from Fishbase
presence <- readRDS("data/Presence_data_Fishbase_Etelis_coruscans.RDS")

## Keep only presence data
presence <- presence[presence$Etelis_coruscans == 1, ]
write.csv(presence, "data/Presence_data_Fishbase_Etelis_coruscans.csv", row.names = F, quote = F)

## Transform as sf object
presence_sf <- st_as_sf(presence, 
                       coords = c("Longitude", "Latitude"), 
                       crs = "WGS84", remove = FALSE)

## Function to shift points lat < 0 to the other side of the map
shift = function(x) {
  geom = st_geometry(x)
  st_geometry(x) = st_sfc(
    lapply(seq_along(geom), function(i) {
      geom[[i]][1] = ifelse(geom[[i]][1] < 0, geom[[i]][1] + 360, geom[[i]][1])
      return(geom[[i]])
    })
    , crs = st_crs(geom)
  )
  return(x)
}

## Mapview plot
mapview(presence_sf)
mapview(shift(presence_sf))

## GGplot map
## Charge global country polygons
world <- ne_countries(scale = "medium", returnclass = "sf")
# worldshift <- st_shift_longitude(world)

ggplot(data = world) +
  geom_sf() +
  geom_point(data = presence_sf, aes(x = Longitude, y = Latitude), color="red", size=1) +
  theme_void()
# coord_sf(xlim = c(-88, -78), ylim = c(24.5, 33), expand = FALSE)



# -----------  Center map around Pacific ----------------------------------------
######## NOT YET FINISHED ######## 

mp1 <- fortify(map(fill=TRUE, plot=FALSE))
mp2 <- mp1
mp2$long <- mp2$long + 360
mp2$group <- mp2$group + max(mp2$group) + 1
mp <- rbind(mp1, mp2)

ggplot() + 
  geom_path(aes(x = long, y = lat, group = group), data = mp) + 
  geom_point(data = shift(presence_sf), aes(x = Longitude, y = Latitude), color="red", size=1) +
  scale_x_continuous(limits = c(-100, 260))




## ---- map sampling sites ----

# get bathymetry data
Bathy <- 
  getNOAA.bathy(lon1 = -180, lon2 = 180,
                lat1 = -30, lat2 = 30,
                resolution = 10,
                antimeridian = TRUE, # to be centered around pacific
                keep = TRUE)

temp <- as.data.frame(table(data_samples$site))
colnames(temp) <- c("site", "N")

coord_site2 <-
  coord_site %>% 
  rownames_to_column("site") %>% 
  left_join(temp)
# write_csv(coord_site2, "coord_sites_N.csv")


