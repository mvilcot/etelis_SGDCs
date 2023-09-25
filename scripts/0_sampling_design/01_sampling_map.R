# ---- setup sites data ----

## get bathymetry data
Bathy <- 
  getNOAA.bathy(lon1 = -180, lon2 = 180,
                lat1 = -30, lat2 = 30,
                resolution = 10,
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


# ---- plot map ----

## transform to sf ----
coord_site_sf <- st_as_sf(coord_siteN, 
                        coords = c("longitude", "latitude"), 
                        crs = "WGS84", remove = FALSE)

data_sites_sf <- st_as_sf(data_sites, 
                          coords = c("Longitude_approx", "Latitude_approx"), 
                          crs = "WGS84", remove = FALSE)

## mapview ----
mapview(shift(data_sites_sf), zcol = "site", cex = 1) +
  mapview(shift(coord_site_sf), zcol = "site", cex = "N", label = "N") 

## charge countries ----
world <- ne_countries(scale = "medium", returnclass = "sf")

# ## !!! ggplot ----
# worldshift <- st_shift_longitude(world)
# 
# ggplot(data = worldshift) +
#   geom_sf() +
#   geom_point(data = shift(data_sites_sf), aes(x = Longitude_approx, y = Latitude_approx), color="black", size=1) +
#   geom_point(data = shift(coord_site_sf), aes(x = longitude, y = latitude), color="red", size=1)
#   # scale_x_continuous(limits = c(0, 360))



