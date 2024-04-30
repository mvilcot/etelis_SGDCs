# Wrapper script for the Etelis coruscans analyses

# ---- load libraries ----
## global
library(tidyverse)      # the beautiful beautiful tidyverse
library(ecodist)        # MRM
library(glue)           # use braces for PXA axis

## plot
library(ggplot2)        # plots
library(ggh4x)          # nested facets
library(patchwork)      # multiple plot by hand
library(gridExtra)      # multiple plot from list
library(viridis)        # palette
library(ggrepel)        # text labels on plots

## spatial
library(terra)          # raster-type package
library(tidyterra)      # terra ggplot
library(sdmpredictors)  # bio-oracle portal
library(geodist)        # geographic distance
library(marmap)         # bathymetry data
library(gdistance)
# library(sf) 
## :::::: --> sf causes an issue on distance matrix when loaded !!! #####
# check in which script it is necessary....

## taxonomy
library(fishtree)       # fish tree of life 
library(betapart)       # jaccard diversity
library(ape)

## genetics
library(dartR)          # amazing wrapper for pop genetics analysis
library(adegenet)
library(hierfstat)
library(mmod)
library(pegas)          # genind2loci
library(jacpop)         # generate_pw_jaccard
library(LEA)            # SNMF

## phylogenetics
library(picante)        # pd, mpd
library(phytools)       # midpoint.root
library(ape)



# ---- load data ----

## genetic data ----
data_samples <- read_csv("data/metadata_samples.csv")

for (level in c("site", "station")){ # order levels
  data_samples[[level]] <- 
    data_samples[[level]] %>% 
    ordered(levels = unique(data_samples[order(data_samples$order),][[level]])) %>% 
    droplevels()
}


## presence data ----
data_Etelis <- readRDS("data/Presence_data_Fishbase_Etelis_coruscans.RDS")


## taxonomy data ----
data_taxo <- read_csv("data/Taxonomy_Fishbase.csv")


## depth data ----
data_depth <- read_csv("data/data_species_depth_range_teleo.csv")


## spatial data ----
data_stations <- read_csv("data/metadata_stations.csv")

for (level in c("site", "station")){ # order levels
  data_stations[[level]] <- 
    data_stations[[level]] %>% 
    ordered(levels = unique(data_stations[order(data_stations$order),][[level]])) %>% 
    droplevels()
}

# get mean coordinates by site
data_sites <- 
  data_stations %>% 
  group_by(site) %>% 
  summarise(longitude=mean(Longitude_approx),
            latitude=mean(Latitude_approx),
            number_samples=sum(number_samples)) %>% 
  arrange(factor(site, levels(data_samples$site)))


# ---- personalised plot ----
level = "site"
color_perso <- c(viridis(length(levels(data_samples[[level]]))))
names(color_perso) <- levels(data_samples[[level]])
LABELS <- gsub('_', ' ', names(color_perso))
LABELS <- gsub('W Australia', 'Western Australia', LABELS)
  


# ---- functions ----
## Melt distance matrix
melt.dist <- function(distmat, metric) {
  if(class(distmat)[1] == "dist") {distmat <- as.matrix(distmat)}
  distmat[upper.tri(distmat, diag = T)] <- NA
  distmat <- 
    as.data.frame(distmat) %>% 
    rownames_to_column(paste0(level, "1")) %>% 
    pivot_longer(cols = -paste0(level, "1"), 
                 names_to = paste0(level, "2"), 
                 values_to = metric) %>% 
    na.omit()
  
  return(distmat)
}


## Read a genlight, with personalized parameters
read.genlight <- function(filters = "missind_callrate0.70_maf0.05",
                          level = "site",
                          site2drop = NULL,
                          site2keep = NULL,
                          station2drop = NULL,
                          station2keep = NULL){
  
  # read genlight  
  genlight <- readRDS(paste0("intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, ".RDS"))
  
  # set population to site
  genlight@pop <- 
    genlight@other[["ind.metrics"]][["site"]] %>%
    droplevels()
  
  # drop sites
  if (!is.null(site2drop)){
    genlight <- gl.drop.pop(genlight, pop.list = site2drop, recalc = T, mono.rm = T)
  }
  
  # keep populations
  if (!is.null(site2keep)){
    genlight <- gl.keep.pop(genlight, pop.list = site2keep, recalc = T, mono.rm = T)
  }
  
  # set population to station
  genlight@pop <- 
    genlight@other[["ind.metrics"]][["station"]] %>%
    droplevels()
  
  # drop sites
  if (!is.null(station2drop)){
    genlight <- gl.drop.pop(genlight, pop.list = station2drop, recalc = T, mono.rm = T)
  }
  
  # keep populations
  if (!is.null(station2keep)){
    genlight <- gl.keep.pop(genlight, pop.list = station2keep, recalc = T, mono.rm = T)
  }
  
  # set population to final level
  genlight@pop <- 
    genlight@other[["ind.metrics"]][[level]] %>%
    droplevels()

  # return
  return(genlight)  
}



## Shift longitude [-180,180] to [0, 360]
# for both df or sf
shift.lon = function(x) {
  longitudecolumn <- grep("Longitude|longitude", colnames(x))
  x[[longitudecolumn]] <- sapply(x[[longitudecolumn]], function(i){ifelse(i < 0, i + 360, i)})

  if("sf" %in% class(x)){
    geom = st_geometry(x)
    st_geometry(x) = st_sfc(
      lapply(seq_along(geom), function(i) {
        geom[[i]][1] = ifelse(geom[[i]][1] < 0, geom[[i]][1] + 360, geom[[i]][1])
        return(geom[[i]])
      }),
      crs = st_crs(geom)
    )
  }
  return(x)
}




## remove location from distance matrix
mat.subset <- function(distmat, location){
  if("dist" %in% class(distmat)){
    distmat <- as.matrix(distmat)
    distmat <- distmat[!grepl(location, rownames(distmat)),
                       !grepl(location, colnames(distmat))]
    distmat <- as.dist(distmat)
  }
  if("data.frame" %in% class(distmat)){
    distmat <- distmat[!grepl(location, distmat$site),]
  }
  distmat
}






# ---- arborescence ----
dir.create("intermediate/", showWarnings = F)
dir.create("intermediate/0_sampling_design/", showWarnings = F)
dir.create("intermediate/1_genetic_diversity/", showWarnings = F)
dir.create("intermediate/2_species_diversity/", showWarnings = F)
dir.create("intermediate/3_distance_metrics/", showWarnings = F)
dir.create("intermediate/4_continuity/", showWarnings = F)

dir.create("results/", showWarnings = F)
dir.create("results/0_sampling_design/", showWarnings = F)
dir.create("results/1_genetic_diversity/", showWarnings = F)
dir.create("results/2_species_diversity/", showWarnings = F)
dir.create("results/3_distance_metrics/", showWarnings = F)
dir.create("results/4_continuity/", showWarnings = F)



# ---- scripts ----
# 1. ...
# source("scripts/1_genetic_diversity/11_snp_fitering.R")




