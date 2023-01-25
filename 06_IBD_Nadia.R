
################################   Genetic distance (FST) ~ Geographic distance  #############################

library("dartR")

# ----------- Read dataset -----------------------------------------------------

## Set parameters :

# sites <- "allsites"
 sites <- "noCOCO" 
# sites <- "noSEYCH_noCOCO"
 filters <- "missind_callrate0.70_maf0.05"
# filters <- "missind_callrate0.70_maf0.05_pcadapt"

genlight <- readRDS(paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_", sites, ".RDS"))
genlight

## Create genlight at station level :
genlight_stations <- genlight
genlight_stations@pop <- genlight_stations@other[["ind.metrics"]][["station"]] 
# Remove populations with monomorphic loci :
genlight_stations <- gl.drop.pop(genlight_stations, pop.list = names(which(table(genlight_stations@pop) < 2)), recalc = T, mono.rm = T)
genlight_stations <- gl.drop.pop(genlight_stations, pop.list = names(which(table(genlight_stations@pop) < 5)), recalc = T, mono.rm = T)


# ----------- Genetic distance -----------------------------------------------------

## Compute Pairwise Fst between all stations :
FST <- gl.fst.pop(genlight_stations)
View(FST)
fst_table <- FST$Fsts  
View(fst_table)
# keep only the lower part of the matrix (values under the diagonal):
library("gdata")
fst_table <- lowerTriangle(fst_table, diag=FALSE, byrow=FALSE)


# ----------- Geographic distance -----------------------------------------------------

# Subset with only longitude latitude :
coordinates <- data.frame(tableCT$Row.names,tableCT$Longitude_approx,tableCT$Latitude_approx) # depends on code from script "05_Distance_to_CT" to have "tableCT"
View(coordinates)
# Rename columns :
colnames(coordinates) <- c("Station","Longitude", "Latitude")


# Compute euclidean distances in degrees:
library("geodist")
dist(coordinates, diag=T, upper=T)
# geodesic distances: 
geodist <- geodist(coordinates,  measure = "geodesic")
View(geodist)
# Rename titles for rows and columns:
rownames(geodist) <- coordinates$Station
colnames(geodist) <- coordinates$Station
library("gdata")
geodist <- lowerTriangle(geodist, diag=FALSE, byrow=FALSE)
plot(fst_table ~ geodist, ylab = "Pairwise FST", xlab = "Geographical distance (m)") # => try without seychelles


        # ALTERNATIVE TO COMPUTE DISTANCES:
        # Convert dataframe to an sf object by using st_as_sf() and specifying the column names of the coordinates and the coordinate reference system:
        library(sf)
        coordinates_sf <- st_as_sf(coordinates, coords = c("Longitude", "Latitude"), crs = 4326)
        # Create a matrix of the distances between all the points :
        st_distance(coordinates_sf)
        dist_matrix <- st_distance(coordinates_sf)
        View(dist_matrix)
        # Rename titles for rows and columns:
        rownames(dist_matrix) <- coordinates$Station
        colnames(dist_matrix) <- coordinates$Station

        
## Without Seychelles:
# /!\ FIRSTLY reload the genlight with filter: sites <- "noSEYCH_noCOCO" and run again everything.
# Subset with only longitude latitude:
coordinates <- data.frame(tableCT$Row.names,tableCT$Longitude_approx,tableCT$Latitude_approx)
# Rename columns:
colnames(coordinates) <- c("Station","Longitude", "Latitude")
# Delete Seychelles:
coord_noSEYCH <- coordinates[-19,]
View(coord_noSEYCH)
# Compute euclidean distances in degrees:
library("geodist")
dist(coord_noSEYCH, diag=T, upper=T)
# geodesic distances: 
geodist_noSEYCH <- geodist(coord_noSEYCH,  measure = "geodesic")
library("gdata")
geodist_noSEYCH <- lowerTriangle(geodist_noSEYCH, diag=FALSE, byrow=FALSE)
plot(fst_table ~ geodist_noSEYCH, ylab = "Pairwise FST", xlab = "Geographical distance (m)") # without seychelles



