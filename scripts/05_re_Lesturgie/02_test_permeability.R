# source("scripts/05_re_Lesturgie/01_crea_raster.R")
library(ade4)


## ---- setup data ----

level = "site"
gd_beta <- readRDS(paste0("results/01_genetic_diversity/gd_list_pairwise_", level, ".RDS"))

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
  
  
  
  fst <- mat_fst/(1-mat_fst)
  
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



## ---- apply ----
coords <- 
  coord_sites %>% 
  dplyr::select(c(longitude, latitude)) %>% 
  as.matrix()

perm <- cell.resistance(raster = raster_final_crop, 
                        param_file = param_file,
                        mat_fst = gd_beta$Fst,
                        coords = coords,
                        IBR = F,
                        LC = T
                        )

perm$best_permeability

plot(perm$sss[,1], perm$sss[,3])
plot(perm$sss2[,2], perm$sss2[,3])


## ---- get least cost distance mat from resistance ---
v <- as.numeric(raster$layer@data@values)
q <- as.numeric(v)
q[which(v==500)] <- 1 # ocean
q[which(v==1000)] <- 23 # seamounts
q[which(v==2000)] <- 100 # coral

raster_new <- setValues(raster_final_crop, q)
tr1 <- transition(raster_new, mean, directions=8)  
tr2 <- geoCorrection(tr1, type="r")
mat_IBR <- costDistance(tr2, coords)
mantel.randtest(mat_IBR, gd_beta$Fst)

write.csv(as.matrix(mat_IBR), "intermediate/05_re_Lesturgie/least_cost_distance_IBR_23_100.csv")


