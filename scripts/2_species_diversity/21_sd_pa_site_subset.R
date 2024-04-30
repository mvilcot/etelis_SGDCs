
# This script creates a PA matrix for each community subset for each station
# The output is a list of PA matrices per community delineation across stations

# ---- load ----
data_PA <- readRDS("data/PA_Mat_GaspObis.RDS")


# ---- define stations ----
# extract sample sites from meta data
vect_stations <- 
  data_stations %>% 
  dplyr::select(c(station, Longitude_approx, Latitude_approx)) %>%
  vect(geom = c("Longitude_approx", "Latitude_approx"), crs = "WGS84")
writeVector(vect_stations, "intermediate/0_sampling_design/vect_stations.gpkg")


# add 100km buffer around each sample site
vect_stations <- 
  vect_stations %>% 
  buffer(width = 50000)


# split PA table into chuncks of 250 species (terra rast too heavy otherwise)
PA_split <- split.default(data_PA[,-c(1,2)], ceiling(seq_along(data_PA[,-c(1,2)]) / 250))

# add lon and lat info
PA_split <- lapply(PA_split, function(table){cbind(data_PA[,c(1,2)], table)})


# function to convert PA lat/lon to PA by location
pa.lonlat2site <- function(PAdf){
  # convert PA to raster and intersect of sample sites and pa raster
  PAdf <- 
    PAdf %>%
    rast(crs = "WGS84") %>% 
    terra::extract(vect_stations,
                   xy = TRUE,
                   ID = FALSE,
                   bind = TRUE,
                   fun=max,
                   na.rm=TRUE)
  
  # PA for each station
  PAdf <-
    as.data.frame(PAdf) %>%
    pivot_longer(cols = -c(station), # pivot longer to handle
                 names_to = "species") %>% 
    group_by(species, station) %>%  # group by species and site
    summarise(pa = max(value, na.rm = TRUE), .groups = "keep") %>%  # count presence within sites 
    pivot_wider(names_from = "species", # reformat to pa mat format
                values_from = "pa")
  
  return(PAdf)
}

# transform each PA into PA by station
PAstation_split <- lapply(PA_split, pa.lonlat2site)

# re-join all
PAstation <- plyr::join_all(dfs = PAstation_split, type = 'left', by = 'station')

# order
PAstation <-
  PAstation %>% 
  arrange(factor(station, levels = levels(data_samples$station)))

PAstation$station <- ordered(PAstation$station)

# save
saveRDS(PAstation, "intermediate/2_species_diversity/PA_Mat_GaspObis_station.RDS")


# ---- PA by site ----

PAsite <-
  data_stations %>% 
  dplyr::select(site, station) %>% 
  left_join(PAstation, by = "station") %>% 
  pivot_longer(cols = -c(site, station), names_to = "species") %>% 
  group_by(site, species) %>% 
  summarise(pa = max(value, na.rm = TRUE), .groups = "keep") %>%  # count presence within all stations of a site 
  pivot_wider(names_from = "species", values_from = "pa")


# # order
# PAsite <-
#   PAsite %>% 
#   arrange(factor(site, levels = levels(data_samples$site)))

saveRDS(PAsite, "intermediate/2_species_diversity/PA_Mat_GaspObis_site.RDS")



