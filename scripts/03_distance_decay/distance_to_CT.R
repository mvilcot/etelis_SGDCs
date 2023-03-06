
## ---- Geographic distances to CT ----

# coordinates of the center of the Coral Triangle
center_CT <- 
  data.frame(longitude=133.679826, 
             latitude=-1.307436)

# compute distance from the different stations to the CT (in meters)
coord_pop$dist_to_CT <- 
  pointDistance(coord_pop[,6:5], center_CT, lonlat=TRUE) 

# Create new dataframe with distance to CT values:
dist_CT <- pointDistance(coord_pop[,6:5], center_CT, lonlat=TRUE) 
# dist_CT
CT <- data.frame(dist_CT) ; print(CT)


## ---- Genetic diversity metrics ----
sites <- "noCOCO" 
filters <- "missind_callrate0.70_maf0.05"
genlight <- readRDS(paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_", sites, ".RDS"))
genlight

# Create genlight at station level :
genlight_stations <- genlight
genlight_stations@pop <- genlight_stations@other[["ind.metrics"]][["station"]] 
# Remove populations with monomorphic loci :
genlight_stations <- gl.drop.pop(genlight_stations, pop.list = names(which(table(genlight_stations@pop) < 2)), recalc = T, mono.rm = T)
# Create dataframe with genetic diversity estimates for each station :
heterozyg <- gl.report.heterozygosity(genlight_stations) 

# View(heterozyg)
# Rename row names of coord_pop dataframe (replace numbers by full name of stations) :
rownames(coord_pop) <- coord_pop$Full.name
# Create new dataframe combining metadata (coord_pop) and diversity estimates (heterozyg) :
test <- merge(coord_pop, heterozyg, by = "row.names", all = T)  # all = T pour garder tous sites meme les sites avec 1 seul inidivus -> il mettra NA
# merge()  by = "row.names" pour fusionner 2 tableaux qui ont juste une colonne en commun (mm nom)
# View(test)




genind <- gl2gi(genlight_stations)
Hs <- Hs(genind)
test2 <- merge(as.data.frame(Hs), coord_pop, by = "row.names")

plot(test2$Hs ~ test2$dist_to_CT, xlab = "Distance to Coral Triangle (m)", ylab = "FIS")
text(test2$Hs ~ test2$dist_to_CT, labels = test$Station)





## ---- Genetic diversity metrics ~ Distance to CT ----

plot(test$FIS ~ test$dist_to_CT, xlab = "Distance to Coral Triangle (m)", ylab = "FIS")
text(test$FIS ~ test$dist_to_CT, labels = test$Station)

plot(test$Hs ~ test$dist_to_CT , xlab = "Distance to Coral Triangle (m)", ylab = "He")
text(test$He ~ test$dist_to_CT, labels = test$Station)

barplot(test$Ho ~ test$Station, xlab = "Station", ylab = "Ho")

plot(test$number_samples ~ test$He) # so try to eliminate sites with less than 5 indiv :

genlight_stations <- gl.drop.pop(genlight_stations, pop.list = names(which(table(genlight_stations@pop) < 5)), recalc = T, mono.rm = T)
heterozyg <- gl.report.heterozygosity(genlight_stations) 
rownames(coord_pop) <- coord_pop$Full.name
test <- merge(coord_pop, heterozyg, by = "row.names", all = F) 
View(test)
tableCT <- test[-c(4,11,13,17,21,22,24,25,27,31,32,33),]  # suppr lignes NA, but not necessary..?
View(tableCT)




# test linear models :
model <- lm(test$Ho ~ test$dist_to_CT)
summary(model)
m2 <- lm(test$FIS ~ test$dist_to_CT)
summary(m2)



