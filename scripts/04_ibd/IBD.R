
library("adegenet")
library("hierfstat")
library("geodist")
library("ecodist")
library("vegan")
library("reshape2")

# https://popgen.nescent.org/StartSNP.html
# https://popgen.nescent.org/DifferentiationSNP.html



# ----------- Setup data -------------------------------------------------------

sitesmeta <- read.csv("Data/sites_coordinates.csv", row.names = 3)
sitesmeta$Full.name <- rownames(sitesmeta)

## Read GD and SD distance matrix
SDbeta <- readRDS("Results/SD_beta_sites_betapair.RDS")
SDbeta <- as.matrix(SDbeta$beta.jtu)
GDbeta <- readRDS("Results/GD_beta_sites_noinf2_pairwise_Fst.RDS")
GDbeta <- as.matrix(GDbeta)

# #### !!!! CHANGE SD row names #######
# rownames(SDbeta) <- colnames(SDbeta) <- 
#   sitesmeta[match(rownames(SDbeta), sitesmeta$Pop), ]$Full.name

# 
# ## Keep only sampling sites present in both GD and SD
SDbeta <- SDbeta[rownames(SDbeta) %in% rownames(GDbeta),
                 colnames(SDbeta) %in% colnames(GDbeta)]
GDbeta <- GDbeta[rownames(GDbeta) %in% rownames(SDbeta),
                 colnames(GDbeta) %in% colnames(SDbeta)]

## Order rows alphabetically
SDbeta <- SDbeta[order(rownames(SDbeta)), order(colnames(SDbeta))]
GDbeta <- GDbeta[order(rownames(GDbeta)), order(colnames(GDbeta))]


## Merge two distance matrix into one data.frame
meltSD <- t(combn(colnames(SDbeta), 2))
meltSD <- data.frame(meltSD, jtu=SDbeta[meltSD])

meltGD <- t(combn(colnames(GDbeta), 2))
meltGD <- data.frame(meltGD, Gstpp=GDbeta[meltGD])

betaFULL <- merge(meltSD, meltGD, all = T)
betaFULL$sites <- paste0(betaFULL$X1, "-", betaFULL$X2)
betaFULL <- betaFULL[,c(5,3,4)] #reorder


# ----------- Geographic distance ----------------------------------------------

coord <- sitesmeta[,c(5,4)]
colnames(coord) <- c("longitude", "latitude")
geodist <-  geodist(coord, measure = "geodesic")
colnames(geodist) <- rownames(geodist) <- rownames(coord)
geodist <- geodist[order(rownames(geodist)), order(colnames(geodist))]

## Merge with full distance data frame
meltGEO <- t(combn(colnames(geodist), 2))
meltGEO <- data.frame(meltGEO, geo=geodist[meltGEO])
meltGEO$sites <- paste0(meltGEO$X1, "-", meltGEO$X2)

############ ERREUR ICI ######################
betaFULL <- merge(betaFULL, meltGEO, all = , by = "sites")
betaFULL <- betaFULL[,c(5,3,4)] #reorder



# ----------- Least-cost distance ----------------------------------------------




# ----------- by site --------------------------------------------------------------
lat <- aggregate(sitesmeta$Latitude_approx, list(sitesmeta$Pop), FUN=mean)
colnames(lat) <- c("sites", "latitude")
lon <- aggregate(sitesmeta$Longitude_approx, list(sitesmeta$Pop), FUN=mean) 
colnames(lon) <- c("sites", "longitude")
coordSITE <- merge(lat, lon, by = "sites")
row.names(coordSITE) <- coordSITE$site


geodist <-  geodist(coordSITE, measure = "geodesic")
colnames(geodist) <- rownames(geodist) <- rownames(coordSITE)
geodist <- geodist[order(rownames(geodist)), order(colnames(geodist))]

## Merge with full distance data frame
meltGEO <- t(combn(colnames(geodist), 2))
meltGEO <- data.frame(meltGEO, geo=geodist[meltGEO])
meltGEO$sites <- paste0(meltGEO$X1, "-", meltGEO$X2)

betaFULL <- merge(betaFULL, meltGEO, all = , by = "sites")




# ----------- IBD --------------------------------------------------------------

## MRM
summary(lm(jtu ~ geo, data=betaFULL))
summary(lm(Gstpp ~ geo, data=betaFULL))

## Plot
ggplot(betaFULL) +
  geom_point(aes(geo, jtu))

ggplot(betaFULL) +
  geom_point(aes(geo, Gstpp))
