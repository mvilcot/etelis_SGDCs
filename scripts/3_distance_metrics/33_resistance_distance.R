
dist_mat <-
  readRDS("intermediate/3_distance_metrics/dist_geo_envt.RDS")



#------------------------------------------------------------------------------#
#--------------------------- 1. CREATE RASTER ----------------------------------
#------------------------------------------------------------------------------#

## ---- empty raster ----
raster_empty <- raster(res = 0.6)


## ---- coord of reefs: convert in a raster ----
# reefs coordinates provided by ReefBase (http://www.reefbase.org)
coord_reefs <- read.csv("data/Lesturgie_2023/raster/Reefs/ReefLocations.csv")
coord_reefs <- cbind(coord_reefs$LON, coord_reefs$LAT)
raster_corals <- rasterize(coord_reefs, raster_empty, field=1)

# area_corals <- area(raster_corals)
# area_cells_c <- sqrt(mean(area_corals[]))

# attribute reefs value to 2000
A2 <- raster_corals@data@values
A2[which(is.na(A2)==F)] <- 2000
raster_corals <- setValues(raster_corals, A2)


## ---- coord of seamounts: convert in a raster ----
# downloaded from https://data.unep-wcmc.org/datasets/41
shp_seamounts <- shapefile("data/Lesturgie_2023/raster/Seamounts/Seamounts.shp") 
coord_seamounts <- cbind(shp_seamounts$LONG, shp_seamounts$LAT)
raster_seamounts <- rasterize(coord_seamounts, raster_empty, field=1)

# attribute seamounts value to 1000
A1 <- raster_seamounts@data@values
A1[which(is.na(A1)==F)] <- 1000
raster_seamounts <- setValues(raster_seamounts, A1)


## ---- raster of the world ----
# downloaded at https://www.naturalearthdata.com/downloads/50m-raster-data/50m-natural-earth-1/
raster_world <- raster("data/Lesturgie_2023/raster/NE1_50M_SR_W/NE1_50M_SR_W.tif")

# area_world <- area(raster_world)
# area_cells_w <- sqrt(mean(area_world[]))
# fa <- area_world/area_cells_w

raster_world <- aggregate(raster_world, fact=18, fun=mean)


## ---- merge 3 rasters + pacific centering ----
raster_merge_0 <- merge(raster_corals, raster_seamounts, raster_world)
raster_merge_W <- crop(raster_merge_0, extent(-180, 0, -90, 90))
raster_merge_E <- crop(raster_merge_0, extent(0, 180, -90, 90))   
extent(raster_merge_W) <- c(180, 360, -90, 90)

raster_merge_PA <- merge(raster_merge_W, raster_merge_E)



## ---- attribute values to ocean and land ----
v <- raster_merge_PA$layer@data@values
w <- v
for (i in 1:length(v)){
  if (v[i]<1000){   # if neither coral or seamounts
    if (v[i]<=140){w[i] <- 500} # 500 for ocean (<140m)
    if (v[i]>140){w[i] <- 0.000001} # 0.000001 for land 
  }
}

raster_final <- setValues(raster_merge_PA, w)




## ---- plot reefs + seamounts ----

df_final <- as.data.frame(raster_final, xy = TRUE)
df_final$category <- df_final$layer
df_final$category[df_final$category==0.000001] <- "land"
df_final$category[df_final$category==500] <- "ocean"
df_final$category[df_final$category==1000] <- "seamount"
df_final$category[df_final$category==2000] <- "reef"

ggplot() +
  geom_raster(data = df_final, aes(x = x, y = y, fill = category)) +
  scale_fill_manual(values = c("#233143","#E0EAEB","#EDAE49","#d1495b")) +
  geom_point(data = shift.lon(coord_site), aes(longitude, latitude), color = "darkgreen") +
  theme_void() +
  coord_fixed()
ggsave("results/3_distance_metrics/map_seamounts_corals_ocean.pdf",
       height = 8, width = 16)


## ---- crop to have only IndoP ----
raster_final_crop <- crop(raster_final, extent(35, 280, -35, 35))

cat("Construction raster complete \n")





#------------------------------------------------------------------------------#
#---------------------- 2. BEST PERMEABILITY VALUES ----------------------------
#------------------------------------------------------------------------------#

library(ade4)

## ---- setup data ----

level = "site"
gd_alpha <- read.csv("results/1_genetic_diversity/gd_table_site.csv")
gd_beta <- readRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))


## ---- subset sites ----
loc = "Seychelles"

# beta
for(i in 1:length(gd_beta)){
  gd_beta[[i]] <- as.matrix(gd_beta[[i]])
  gd_beta[[i]] <- gd_beta[[i]][!grepl(loc, rownames(gd_beta[[i]])),
                                 !grepl(loc, colnames(gd_beta[[i]]))]
  gd_beta[[i]] <- as.dist(gd_beta[[i]])
}

# alpha
gd_alpha <- gd_alpha[!grepl(loc, gd_alpha$site),]

# coord
coord_site <- coord_site[!grepl(loc, rownames(coord_site)),]



v <- raster_final_crop$layer@data@values
# hist(v)

# cell values in the raster : 
# 2000 (coral)
# 1000 (seamount)
# 500 (open ocean)
# 1e-06 (land)

param_file <- data.frame(
  c(1,100,1,2000), # coral reef value is investigated from 1 to 100 by 0.01
  c(1,'d',1,1000), # seamount value is investigated from 1 to the coral reef value "d" by 0.01
  c(1,'fixed','fixed',500), # ocean is fixed to 1
  c(1e-6,'fixed','fixed',1e-06)) # land is fixed to 0. 
colnames(param_file) = c('corals','seamounts','ocean','land')
rownames(param_file) = c('value','supp_bound_range','iter_by','raster_cell_value')
param_file



## ---- function ----
cell.resistance <- function(raster, param_file, mat_fst, coords, IBR=F, LC=T){
  h <- c()
  d <- c()
  v <- as.numeric(raster$layer@data@values)
  
  # 
  for (col in 1:ncol(param_file)){
    if (is.numeric(param_file[2,col])) {h <- c(h,col)} # h for set upper values
    else if (param_file[2,col]=='d') {d <- c(d,col)} # d for dependant values
    else {
      v[which(v==param_file[4,col])] <- param_file[1,col] # for fixed values, replace raster value by resistance value
    }
  }
  
  lh <- length(h) # how many variable needs to vary
  
  byl=1
  val <- list()
  for(i in 1:length(h)){
    val[[i]] <- seq(from=param_file[1,i],to=param_file[2,i],by=param_file[3,i]) # values to increment
    byl <- (length(val[[i]]))*byl
  }
  sss <- matrix(NA,ncol=(lh+1+length(d)),byl)
  if(lh==2){
    jj <- c()
    for(i in 1:length(val[[1]])){
      jj <- c(jj,rep(val[[1]][i],length(val[[2]])))
    }
    jj2 <- c(rep(val[[2]],length(val[[1]])))
    sss[,1] <- jj
    sss[,2] <- jj2
  }else{
    sss[,1]=val[[1]]
    sss[,2]=c(rep(as.numeric(as.character(param_file[1,d])),length(val[[1]])))
  }
  
  
  
  # fst <- mat_fst/(1-mat_fst)
  fst <- mat_fst
  
  # list of positions of each variable (reefs, seamounts)
  vlist <- list()
  for(i in 1:lh){
    vlist[[i]] <- which(v==param_file[4,h[i]])
  }
  if(length(d)>0){
    vlist[[lh+1]] <- which(v==param_file[4,d])
  }
  
  # iteration for each couple of values (k)
  q <- as.numeric(v)
  for(k in 1:nrow(sss)){
    q[vlist[[1]]] <- sss[k,1] # set to initial value
    q[vlist[[2]]] <- sss[k,2] # set to initial value
    raster_new <- setValues(raster, q)
    tr1 <- transition(raster_new, mean, directions=8)  
    tr2 <- geoCorrection(tr1, type="r")
    if(IBR==T){
      temp <- commuteDistance(tr2,coords)
    }else if(LC==T){
      temp <- costDistance(tr2, coords)
    }
    mantel <- mantel.randtest(temp,fst)
    sss[k,3] <- mantel$obs
    #print(sss[k,])
    if(length(d)>0){
      cat('Variable 1: iteration ',k, "/", nrow(sss), "\n")
    }else{cat('Iteration ',k, "/", nrow(sss), "\n")}
    
  }
  best_permeability <- sss[which(sss[,3]==max(sss[,3])),]
  best_permeability <- t(data.frame(best_permeability))
  
  if(length(d)>0){
    byl2 <- ((best_permeability[1,1]/as.numeric(as.character(param_file[1,d])))/as.numeric(as.character(param_file[3,d])))
    
    seq2 <- seq(as.numeric(as.character(param_file[1,d])),best_permeability[1,1],as.numeric(as.character(param_file[3,d])))
    byl2 <- length(seq2)
    sss2 <- matrix(NA,ncol=(lh+1+length(d)),byl2)
    sss2[,1] <- rep(best_permeability[1,1],byl2)
    sss2[,2] <- seq2
    for(k in 1:nrow(sss2)){
      q[vlist[[1]]] <- sss2[k,1]
      q[vlist[[2]]] <- sss2[k,2]
      raster_new <- setValues(raster, q)
      tr1 <- transition(raster_new, mean, directions=8)  
      tr2 <- geoCorrection(tr1, type="r")
      if(IBR==T){
        temp <- commuteDistance(tr2,coords)
      }else if(LC==T){
        temp <- costDistance(tr2, coords)
      }
      
      mantel <- mantel.randtest(temp,fst)
      sss2[k,3] <- mantel$obs
      #print(sss[k,])
      cat('Variable 2: iteration ',k, "/", nrow(sss2), "\n")
    }
    best_permeability <- sss2[which(sss2[,3]==max(sss2[,3])),]
    best_permeability <- t(data.frame(best_permeability)) 
  }
  
  if(length(d)>0){
    colnames(best_permeability) <- c(paste0("p_",colnames(param_file)[h]),paste0("p_",colnames(param_file)[d]),"R")
  }else{
    colnames(best_permeability) <- c(paste0("p_",colnames(param_file)[h[1]]),paste0("p_",colnames(param_file)[h[2]]),"R")
  }
  return(list(best_permeability = best_permeability, 
              sss = sss,
              sss2 = sss2))
  
}



## ---- find best permeability value ----

perm <- cell.resistance(raster = raster_final_crop, 
                        param_file = param_file,
                        mat_fst = gd_beta$Fst,
                        coords = as.matrix(coord_site),
                        IBR = F,
                        LC = T
)

perm$best_permeability

plot(perm$sss[,1], perm$sss[,3])
plot(perm$sss2[,2], perm$sss2[,3])





## ---- 3. LEAST COST DISTANCE ----
v <- as.numeric(raster_final_crop$layer@data@values)
q <- as.numeric(v)
q[which(v==500)] <- 1 # ocean
q[which(v==1000)] <- as.data.frame(perm$best_permeability)$p_seamounts # seamounts
q[which(v==2000)] <- as.data.frame(perm$best_permeability)$p_corals # coral

raster_new <- setValues(raster_final_crop, q)
tr1 <- transition(raster_new, mean, directions=8)  
tr2 <- geoCorrection(tr1, type="r")

# get mean coordinates by location
coord_site <- 
  data_sites %>% 
  group_by(site) %>% 
  summarise(longitude=mean(Longitude_approx),
            latitude=mean(Latitude_approx)) %>% 
  arrange(factor(site, levels(data_samples$site))) %>%
  column_to_rownames("site")


dist_mat$leastcost <- 
  gdistance::costDistance(tr2, as.matrix(shift.lon(coord_site)))

mantel.randtest(mat_IBR, gd_beta$Fst)



## ---- export ----

dist_mat %>% 
  saveRDS("intermediate/3_distance_metrics/dist_geo_envt_res17-4.RDS")


