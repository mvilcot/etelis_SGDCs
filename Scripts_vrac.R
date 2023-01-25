

# ----------- Read dataset -----------------------------------------------------

library("adegenet")
library("hierfstat")
library("dartR")
library("mmod")

## Set parameters
sites <- "allsites" 
filters <- "missind_callrate0.70_maf0.05_pcadapt"

genlight <- readRDS(paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_", sites, ".RDS"))
genlight



# ----------- Convert to other formats -----------------------------------------
## Genind format
genind <- dartR::gl2gi(genlight) #same than adeget::df2genind
genind

## Genpop format
genepop <- adegenet::genind2genpop(genind) 
genepop

## Hierfstat format
hierfstatDF <- hierfstat::genind2hierfstat(genind, pop = NULL)



# ----------- Genetic diversity metrics ----------------------------------------

# https://popgen.nescent.org/StartSNP.html
# https://popgen.nescent.org/DifferentiationSNP.html

## Overall statistics
BS <- hierfstat::basic.stats(hierfstatDF, diploid = TRUE, digits = 4)
BS$overall # Statistics over all samples

## G"st
Gst_Hedrick(genind)

## Heterozigosity within population
Hs(genepop, pop = genepop$other$ind.metrics$pop)

## Pairwise G"st distance
pairwise_Gst_Hedrick(genind)
#ou
genet.dist(hierfstatDF, diploid=TRUE, method="Nei87")



#### -------- Least cost distance -----------------------------------------------------------------------------------
# Geographic distance is computed as least-cost path distances, 
# i.e. the length of the shortest path within water between two sites.
# Je t'ai mis un script que j'avais utilisé sur un autre jeu de données, à adapter au Vivandeau donc

library("marmap")

## Read data
coord_pop <- read.csv("Data/Sites_coordinates.csv", row.names = 1) # Population mean coordinates

## Get bathymetry data for the zone
Bathy <- getNOAA.bathy(lon1 = min(coord_pop[,"Longitude"])-5, lon2 = max(coord_pop[,"Longitude"])+5,
                       lat1 = min(coord_pop[,"Latitude"])-5, lat2 = max(coord_pop[,"Latitude"])+5,
                       resolution = 1)
blues <- colorRampPalette(c("red", "purple", "blue", "cadetblue1", "white"))

## Plot 
pdf(file = "Results/04_Bathymetry_plot.pdf", height=7, width=9.5)
plot(Bathy, image = TRUE, bpal = blues(50))
scaleBathy(Bathy, deg = 5, x = "topleft", inset = 5)
points(coord_pop[,"Longitude"], coord_pop[,"Latitude"],
       pch = 21, col = "black", bg = "red", cex = 1.3)
dev.off()

## Compute pairwise within-water distance
# Voir point 2.1 :
# https://cran.r-project.org/web/packages/marmap/vignettes/marmap-DataAnalysis.pdf
Trans <- trans.mat(Bathy, min.depth=-10, max.depth=NULL)
dist_geo <- lc.dist(Trans, coord_pop, res=c("dist","path"))
write.csv(as.matrix(dist_geo), file = "Intermediate/Least_cost_distance.csv")







#### -------- creer une grille -------------------------------------------------
grd <- SpatialGrid(GridTopology(cellcentre.offset=c(-180,-90),cellsize=c(1,1),cells.dim=c(361,181)))

#Creation de la liste des coordonnees de la grille PA
grd_coord <- as(grd,"SpatialPolygons")

a <- list()
for(i in 1:length(grd_coord@polygons)){
  a[[i]] <- grd_coord@polygons[[i]]@labpt
} # end of i

coords <- do.call(rbind,a)
grille <- data.frame(ID=seq(1,dim(coords)[1],1),X_ENTIER=coords[,1],Y_ENTIER=coords[,2])

#Switch from grid to raster
grd <- raster(grd); grd[] <-0

# Affiliation d'une projection
grd@crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Inventaire des fichiers shp dans le repertoire de stockage des polygones
list_species <- do.call(rbind,strsplit(list.files(getwd(),pattern=".shp"),".shp"))[,1]

cat("### Making the Matrice P/A ###","\n")

res <- mclapply(1:length(list_species),mc.cores=mc.cores,function(x){
  cat ("i=",x,"\n")
  get_vect_species(sp=list_species[x],rast_grid=grd)
})

mat_PA <- cbind(grille[,2:3],do.call(cbind,res))
colnames(mat_PA)[-c(1:2)] <- list_species


#### -------- déterminer les quantiles pour la répartition bathymétrique -------

for(file_name in liste_teleo){
  data_sp <- read.csv(file_name,row.names=1,dec=".")
  #if (nrow(data_sp>=10)){
  coord_test <- paste(round(data_sp$decimalLongitude,5),round(data_sp$decimalLatitude,5),sep="_")
  if(length(unique(coord_test))<=20){
    res_bath <- c(NA,NA)
  } else{ sp_ext <- extent( cbind(data_sp$decimalLongitude,data_sp$decimalLatitude))
  r <- crop(r_bathy,sp_ext)
  # Sampling point protocol
  cells <-cellFromXY(r, cbind(data_sp$decimalLongitude,data_sp$decimalLatitude))
  data_cells <- r[][cells]
  bathy <- data_cells[which(data_cells<=0)]
  number_occurrence<- length(bathy)
  rm(r)
  # Bathy min and max
  #I get the 1st quartile and the 3rd quartile of the bathymetric distribution
  quanti<- quantile(bathy, probs = seq(0,1,0.01))
  res_bath= c(file_name, bathy_min=abs(quanti[100]),bathy_max=abs(quanti[2]), number_occurrence)
  names(res_bath) <- c("species", "bathy_min","bathy_max", "number_occurrence")
  sp_depth_teleo<- rbind(sp_depth_teleo, res_bath)
  }
}







