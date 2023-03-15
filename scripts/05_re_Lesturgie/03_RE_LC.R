#--------------------------------------------------#
#-------------- >>> LEAST COST <<< ----------------#
#--------------------------------------------------#

## ---- PARAMETERS FOR PERMEABILITY ---- 
value_ocean <- 1
value_coral <- 100
value_seamount <- 23

## ---- CREATING RASTER ---- 
source("scripts/05_re_Lesturgie/01_crea_raster.R")
coord_cells <- coordinates(raster_final_crop)


## ---- ATTRIBUTING PERMEABILITY VALUES ---- 

#500 (ocean)
#2000 (coraux)
#1000 (Seamounts)

mappa <- raster_final_crop

x <- mappa$layer@data@values
y <- x
for (i in 1:length(x)){
  if (x[i]==0.000001){}
  if (x[i]==500){y[i] <- value_ocean}
  if (x[i]==1000){y[i] <- value_seamount}
  if (x[i]==2000){y[i] <- value_coral}
}

mappaLC1 <- setValues(mappa, y)


## ---- CREATING TRANSITION MATRIX ---- 
translay <- transition(mappaLC1, mean, directions = 8)  
translay <- geoCorrection(translay, type = "r") 


## ---- LEAST COST CORRELATION WITH GENETIC DIVERSITY ---- 
mappa_corr_least_cost <- mappaLC1

n_pop <- length(coord_sites[,1])
cell_tot <- length(mappaLC1@data@values)
for (cell in 1:cell_tot) {
  if (mappaLC1@data@values[cell]==0.000001 ) {mappa_corr_least_cost@data@values[cell] <- NA}
  else {
    vett_dist <- matrix(rep(NA,n_pop), ncol=1, nrow=n_pop)
    for (j in 1:n_pop) {
      vett_dist[j,1]  <-  costDistance(translay, as.matrix(rbind(coord_cells[cell,], coord_sites[j,1:2])))
    }
    mappa_corr_least_cost@data@values[cell] <- cor(coord_sites[,3],vett_dist)
  }
  if (cell %% 100 == 0) {
    print(paste0("Computing Least Cost Distance correlations: ",cell,"/",cell_tot," ...cells"))
  }
}

saveRDS(mappa_corr_least_cost, "intermediate/05_re_Lesturgie/mappa_corr_least_cost.RDS")

pdf(paste0("results/05_re_Lesturgie/RE_map_LC_IP_S3_", value_coral, "_", value_seamount, ".pdf"))
plot(mappa_corr_least_cost)
points(coord_sites[,1],coord_sites[,2], pch=20, cex=.7)
dev.off()



mappa_corr_least_cost2 <- disaggregate(mappa_corr_least_cost, 18)

pdf(paste0("results/05_re_Lesturgie/RE_map_LC_IP_S3_qx18", value_coral, "_", value_seamount, ".pdf"))
plot(mappa_corr_least_cost2)
points(coord_sites[,1],coord_sites[,2], pch=20, cex=.7)
dev.off()



save.image(paste0("intermediate/05_re_Lesturgie/RE", value_coral, "_", value_seamount, ".RData"))


