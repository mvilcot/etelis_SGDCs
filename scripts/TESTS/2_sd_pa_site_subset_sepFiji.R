
# This script creates a PA matrix for each community subset for each station
# The output is a list of PA matrices per community delineation across stations

## ---- load ----
data_PA <- readRDS("data/PA_Mat_GaspObis.RDS")


# ## ---- define stations ----
# # extract sample sites from meta data
# sites_vect <- 
#   data_sites %>% 
#   dplyr::select(c(station, Longitude_approx, Latitude_approx)) %>%
#   vect(geom = c("Longitude_approx", "Latitude_approx"), crs = "WGS84")
# # writeVector(sites_vect, "intermediate/00_sampling_sites/sites_vect.gpkg")
# 
# 
# # add 100km buffer around each sample site
# sites_vect <- 
#   sites_vect %>% 
#   buffer(width = 50000)
# # writeVector(sites_vect, "intermediate/00_sampling_sites/sites_vect_buffer.gpkg")
# 
# 
# # create raster of PA grid for Qgis 
# PA_sub <- 
#   data_PA[,1:2] %>% 
#   cbind(1) %>% 
#   rast(crs = "WGS84")
# # writeRaster(PA_sub, "intermediate/00_sampling_sites/matrix_PA.gpkg")
# 
# 
# 
# ## ---- PA by station ----
# 
# # # reduce PA data to species of interest (too heavy otherwise)
# # list_commFULL <-
# #   list_communities %>%
# #   reduce(append) %>%
# #   unique()
# 
# # split PA table into chuncks of 250 species (terra rast too heavy otherwise)
# PA_split <- split.default(data_PA[,-c(1,2)], ceiling(seq_along(data_PA[,-c(1,2)]) / 250))
# 
# # add lon and lat info
# PA_split <- lapply(PA_split, function(table){cbind(data_PA[,c(1,2)], table)})
# 
# 
# # function to convert PA lat/lon to PA by site
# pa.lonlat2site <- function(PAdf){
#   # convert PA to raster and intersect of sample sites and pa raster
#   PAdf <- 
#     PAdf %>%
#     rast(crs = "WGS84") %>% 
#     terra::extract(sites_vect,
#                    xy = TRUE,
#                    ID = FALSE,
#                    bind = TRUE,
#                    fun=max,
#                    na.rm=TRUE)
#   
#   # PA for each station
#   PAdf <-
#     as.data.frame(PAdf) %>%
#     pivot_longer(cols = -c(station), # pivot longer to handle
#                  names_to = "species") %>% 
#     group_by(species, station) %>%  # group by species and site
#     summarise(pa = max(value, na.rm = TRUE), .groups = "keep") %>%  # count presence within sites 
#     pivot_wider(names_from = "species", # reformat to pa mat format
#                 values_from = "pa")
#   
#   return(PAdf)
# }
# 
# # transform each PA into PA by station
# PAstation_split <- lapply(PA_split, pa.lonlat2site)
# 
# # re-join all
# PAstation <- plyr::join_all(dfs = PAstation_split, type = 'left', by = 'station')
# 
# # save
# saveRDS(PAstation, "intermediate/02_species_diversity/PA_Mat_GaspObis_station.RDS")



# ## ---- Method 2: Loop on each species (too long) ----
# list_species <- colnames(data_PA)[!colnames(data_PA) %in% c("Longitude", "Latitude")]
# 
# PAstation <- 
#   data.frame(station = sites_vect$station) %>% 
#   arrange(station)
# 
# 
# for (i in 1:length(list_species)){
#   sp <- list_species[i]
#   if(i %% 500 == 0){ print(i)}
#   
#   # subset PA to 1 species
#   data_PAsub <-
#     data_PA %>%
#     select(all_of(c("Longitude", "Latitude", sp)))
#   
#   
#   # extract 1 PA data frame per community by site
#   # subset PA to community, convert to raster
#   # intersect of sample sites and pa raster
#   PA <- 
#     data_PAsub %>%
#     rast(crs = "WGS84") %>% 
#     terra::extract(sites_vect,
#                    xy = TRUE,
#                    ID = FALSE,
#                    bind = TRUE,
#                    fun=max,
#                    na.rm=TRUE)
#   
#   # pa for each station
#   PA <-
#     as.data.frame(PA) %>%
#     pivot_longer(cols = -c(station), # pivot longer to handle
#                  names_to = "species") %>% 
#     group_by(species, station) %>%  # group by species and site
#     summarise(pa = max(value, na.rm = TRUE), .groups = "keep") %>%  # count presence within sites 
#     pivot_wider(names_from = "species", # reformat to pa mat format
#                 values_from = "pa")
#   
#   PAstation <- 
#     PAstation %>% 
#     left_join(PA, by = "station")
#   
# }
# 
# 
# # export
# saveRDS(PAstation, "intermediate/02_species_diversity/PA_Mat_GaspObis_station2.RDS")






# ## ---- Method 3: subset by station (error GASCOYE) ----
# 
# library(sf)
# data_PAsf <- st_as_sf(data_PA, coords = c("Longitude", "Latitude"), crs=4326, remove=FALSE)
# 
# 
# ## Create empty df
# PAstations <- data.frame(matrix(ncol = ncol(data_PA), nrow = 0))
# colnames(PAstations) <- colnames(data_PA)
# 
# # sf::sf_use_s2(FALSE)
# for (i in 1:nrow(data_sites)){
#   site <- data_sites$station[i]
#   lon <- data_sites$Longitude_approx[i]
#   lat <- data_sites$Latitude_approx[i]
#   print(paste(i, site, lon, lat, sep = ", "))
# 
#   PAtemp <- st_crop(data_PAsf,
#                     xmin = round(lon)-0.5, xmax = round(lon)+0.5,
#                     ymin = round(lat)-0.5, ymax = round(lat)+0.5)
# 
#   st_geometry(PAtemp) <- NULL # convert sf to dataframe
#   if(nrow(PAtemp) == 1){
#     PAstations[site, ] <- PAtemp
#   }else{print(paste("!!!!! ERROR !!!!!", site))}
# }
# 
# saveRDS(PAstations, "intermediate/02_species_diversity/PA_Mat_GaspObis_station2.RDS")
# 
# 
# 


## ---- PA by site ----
PAstation <- readRDS("intermediate/02_species_diversity/PA_Mat_GaspObis_station.RDS")

sites_sub <-
  data_sites %>% 
  dplyr::select(site, station)
  
PAsite <-
  sites_sub %>%
  left_join(PAstation, by = "station") %>% 
  pivot_longer(cols = -c(site, station), names_to = "species") %>% 
  group_by(site, species) %>% 
  summarise(pa = max(value, na.rm = TRUE), .groups = "keep") %>%  # count presence within all stations of a site 
  pivot_wider(names_from = "species", values_from = "pa")

saveRDS(PAsite, "intermediate/TEST/PA_Mat_GaspObis_site_sepFiji.RDS")




## ---- PA all stations ----

PAall <-
  PAsite %>% 
  pivot_longer(cols = -c(site), names_to = "species") %>% 
  group_by(species) %>% 
  summarise(pa = max(value, na.rm = TRUE), .groups = "keep") %>%  # count presence within all stations of a site 
  pivot_wider(names_from = "species", values_from = "pa")

saveRDS(PAall, "intermediate/TEST/PA_Mat_GaspObis_allstations_sepFiji.RDS")




