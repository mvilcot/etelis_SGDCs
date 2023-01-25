
library("hierfstat")
library("radiator")
library("dartR")
library("adegenet")
library("mmod")
library("HierDpart")
library("betapart")
library("diveRsity")

## Set parameters
# filters <- "missind_callrate0.70_maf0.05_pcadapt"
filters <- "missind_callrate0.70_maf0.05"

## Read genlight  
genlight <- readRDS(paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_allsites.RDS"))
genlight

## OPTIONAL: Set population to stations 
genlight@pop <- genlight@other[["ind.metrics"]][["station"]]

## Remove populations with less than two individuals
print(names(which(table(genlight@pop) < 2)))
genlight <- gl.drop.pop(genlight, pop.list = names(which(table(genlight@pop) < 2)), recalc = T, mono.rm = T)
genlight

het <- gl.report.heterozygosity(genlight)

testmerge <- merge(sitesmeta, het, by="row.names")




# ----------- Convert to other formats -----------------------------------------
## Genind format
genind <- dartR::gl2gi(genlight) # same than adegenet::df2genind
genind

## Genpop class format
genpop <- adegenet::genind2genpop(genind) 
genpop

## Genepop format (Rousset)
genomic_converter(genlight, output = "genepop")

genepop <- dartR::gl2genepop(genlight, outfile = 'Intermediate/GenepopdartR_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_allsites.genepop_test2.gen', outpath = ".") 
genepop

## Hierfstat format
hierfstatDF <- hierfstat::genind2hierfstat(genind)



# ----------- Genetic diversity metrics ----------------------------------------

## Overall statistics
BS <- hierfstat::basic.stats(hierfstatDF, diploid = TRUE, digits = 4) # why digits=4 ?
BS
BS_o <- BS$overall # Statistics over all samples
# OU :
BS2 <- basic.stats(genind) # same result
BSpop <- apply(BS2$Ho, MARGIN = 2, FUN = mean, na.rm = T)

testmerge$Station <- factor(testmerge$Station, levels=unique(testmerge[order(testmerge$order),]$Station), ordered=TRUE)

ggplot(testmerge) +
  geom_point(aes(Station, Ho)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# wc(genind) # a tester

div <- summary(hierfstatDF)



## Pairwise Fst:
GDbetaFst <- genet.dist(genind, method = "WC84")  # a faire tourner sur ordi plus puissant?
saveRDS(GDbetaGst, "Results/GD_beta_sites_noinf2_pairwise_Fst.RDS")


## G"st
# Gst_Hedrick(genind) # Gst global = 0.045
GDbetaGst <- pairwise_Gst_Hedrick(genind)
saveRDS(GDbetaGst, "Results/GD_beta_stations_noinf2_pairwise_GstppHedrick.RDS")


## Jost Additive framework 
Hs = BS_o[["Hs"]]
Ht = BS_o[["Ht"]]
Hst = (Ht-Hs)/(1-Hs)
Dst <- Ht - Hs


## Jost Multiplicative framework Jst=Jt*Js
Js = 1/(1-Hs)
Jt = 1/(1-Ht)
Jst2 = Jt/Js # or 1/(1-Hst)


##### ICI #####
## TEST PAIRWISE (diveRsity)
test <- diffCalc(infile = "Intermediate/GenepopRadiator_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_sites_noCocos.gen", 
         outfile = "test", fst = T, pairwise = T)

plot(test$pairwise$D ~ test$pairwise$GGst)


test2 <- fastDivPart(infile = "Intermediate/GenepopRadiator_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_sites_noCocos.gen", outfile = "myresults", gp = 3, 
                 pairwise = T, plot = TRUE)




## Jaccard
GDbetaJAC <- HierJd("Intermediate/test.gen", ncode = 3, r = 1, nreg = 1)
saveRDS(GDbetaJAC, "Results/test_HierJd_betaGD_stations_noinf2.RDS")
# !!!!!!!!!!!!!!!!! Error in HierJd("Intermediate/test.gen", ncode = 3) : !!!!!!!!!!!! ####
# argument "r" is missing, with no default
# In addition: There were 27 warnings (use warnings() to see them)

########### ICI ################
## Jaccard Betapart
PAalleles <- as.data.frame(genind@tab)
PAallelesTEST <- PAalleles
PAallelesTEST[PAallelesTEST > 0] <- 1 

PAallelesTEST$site <- genind@pop
PAallelesTESTagg <- aggregate(PAallelesTEST[,-ncol(PAallelesTEST)], by=list(Site=PAallelesTEST$site), FUN=max) ##by site


test <- beta.multi(PAallelesTEST[,-ncol(PAallelesTEST)], index.family="jaccard")


## Hill ????





# ----------- Read Species PA data ------------------------------------------------------------

library("sf")
library("sp")
library("raster")
library("betapart")

## Sites d'echantillonage
metasites <- read.csv("Data/sites_coordinates.csv", row.names = 2)

## Donnees présence-absence OBIS/GASPAR
PA <- readRDS("Data/PA_Mat_GaspObis.RDS")
PAsf <- st_as_sf(PA, coords = c("Longitude", "Latitude"), crs=4326, remove=FALSE)




# ----------- Species presence per station -------------------------------------

SDalpha <- metasites[, 1:4]
SDalpha$SRcommunity <- NA

## Create empty df
PAstations <- data.frame(matrix(ncol = ncol(PA), nrow = 0))
colnames(PAstations) <- colnames(PA)

# sf::sf_use_s2(FALSE)
for (i in 1:nrow(metasites)){
  site <- rownames(metasites)[i]
  
  lon <- metasites$Longitude_approx[i]
  lat <- metasites$Latitude_approx[i]
  print(paste(i, site, lon, lat, sep = ", "))
  
  PAtemp <- st_crop(PAsf,
                    xmin = round(lon)-0.5, xmax = round(lon)+0.5,
                    ymin = round(lat)-0.5, ymax = round(lat)+0.5)
  
  st_geometry(PAtemp) <- NULL # convert sf to dataframe
  if(nrow(PAtemp) == 1){
    PAstations[site, ] <- PAtemp
  }
  else{print(paste("ERROR", site))}
}

write.csv(PAstations, "Intermediate/PA_species_stations.csv", row.names = T, quote = F)
PAstations <- read.csv("Intermediate/PA_species_stations.csv", row.names = 1)



## Combine by site/population

PAstations <- PAstations[,-c(1,2)]
PAsites <- as.data.frame(setNames(replicate(length(PAstations), numeric(0), simplify = F), colnames(PAstations)))

for(site in unique(metasites$Pop)){
  PAtemp <- PAstations[rownames(metasites[metasites$Pop == site,]),]
  PAsites[site,] <- apply(PAtemp, 2, max)
}



# ----------- Species Alpha diversity --------------------------------------------------

SDalpha$SRcommunity <- rowSums(PAstations)





# ----------- Species Beta diversity --------------------------------------------------

## Per station
## ATTENTION? SANS GASCOYE POUR LE MOMENT !!!!!
PAstations2 <- PAstations[, -c(1,2)] 
PAstations2 <- PAstations2[rownames(PAsites2) != "Gascoyne", ] 
SDbeta <- beta.pair(PAstations2, index.family="jaccard")
# saveRDS(SDbeta, "Results/SD_beta_stations_betapair.RDS")


## Per site
PAsites2 <- PAsites[, -c(1,2)] 
PAsites2 <- PAsites2[rownames(PAsites2) != "W_Australia", ] 
SDbeta <- beta.pair(PAsites2, index.family="jaccard")
# saveRDS(SDbeta, "Results/SD_beta_sites_betapair.RDS")



## Visualise
dst <- data.matrix(SDbeta$beta.jtu)
dim <- ncol(dst)

image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, rownames(PAsites2), cex.axis = 0.5, las=3)
axis(2, 1:dim, rownames(PAsites2), cex.axis = 0.5, las=1)

text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", dst), cex=0.6)
