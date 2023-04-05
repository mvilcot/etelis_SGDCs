## ---- setup ----
dist_mat <- list()


## ---- geographic distance ----
# compute distance between locations (in meters)
dist_mat$geodist <- 
  geodist(coord_site, measure = "geodesic") 

colnames(dist_mat$geodist) <- 
  rownames(dist_mat$geodist) <- 
  rownames(coord_site)

dist_mat$geodist <-
  dist_mat$geodist %>% 
  as.dist()

## ---- distance by sea----

# get bathymetry data
Bathy <- 
  getNOAA.bathy(lon1 = -180, lon2 = 180,
                lat1 = -30, lat2 = 30,
                resolution = 10,
                antimeridian = TRUE, # to be centered around pacific
                keep = TRUE)

saveRDS(Bathy, "intermediate/3_distance_metrics/bathymetry.RDS")
color_blues <- 
  colorRampPalette(c("purple", "blue", "cadetblue1", "cadetblue2", "white"))

# save bathymetry map
pdf(file = "results/3_distance_metrics/bathymetry_plot.pdf", height=7, width=14)
plot.bathy(Bathy, 
           step=2000, 
           deepest.isobath = -14000, shallowest.isobath = 0, 
           image = TRUE, bpal = color_blues(20))

scaleBathy(Bathy, 
           deg = 5, x = "topleft", inset = 5)

points(shift.lon(coord_site)$longitude, shift.lon(coord_site)$latitude,
       pch = 21, col = "black", bg = "red", cex = 1.3)

dev.off()

# compute distance by sea
trans <- 
  trans.mat(Bathy, min.depth=-5, max.depth=NULL)

dist_mat$seadist <- 
  lc.dist(trans, shift.lon(coord_site), res = "dist")



## ---- export ----
dist_mat %>% 
  saveRDS("intermediate/3_distance_metrics/dist_geo.RDS")


