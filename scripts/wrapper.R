# Wrapper script for the Etelis coruscans analyses

## ---- load libraries ----
# scripts
library(tidyverse)      # the beautiful beautiful tidyverse
library(pracma)         # findintervals
library(reshape)        # melt
library(ecodist)        # MRM

# plot
library(ggplot2)        # plots
library(wesanderson)    # palette
library(ggh4x)
library(paletteer)      # palette
library(patchwork)      # easy multiple plot
library(gridExtra)      # easy multiple plot

# spatial
library(terra)          # raster-type package
library(sdmpredictors)  # bio-oracle portal
library(geodist)        # geographic distance
library(marmap)         # bathymetry data
library(gdistance)

# taxonomy
library(fishtree)       # fish tree of life 
library(betapart)       # jaccard diversity
library(ape)

# genetics
library(adegenet)
library(dartR)          # SNPs filtering
library(pcadapt)
library(qvalue)
library(OutFLANK)
library(hierfstat)
library(mmod)
library(pegas)          # genind2loci
library(jacpop)         # generate_pw_jaccard

# phylogenetics
library(picante)        # pd, mpd
library(ggtree)         # beautiful tree plot
library(phytools)       # midpoint.root



## ---- load data ----
# genetic data
data_samples <- read.csv("data/metadata_samples.csv")
data_sites <- read.csv("data/metadata_sites.csv")

data_samples <- 
  data_samples %>% 
  left_join(data_sites, by = c("station", "site"))


# species data
data_Etelis <- readRDS("data/Presence_data_Fishbase_Etelis_coruscans.RDS")

# taxonomy data
data_species <- read.csv("data/data_species_depth_range_teleo.csv")
data_species2 <- read.csv("data/data_species.csv")
data_fishtree <- read.csv("data/PFC_taxonomy.csv")
# data_fishbase <- rfishbase::load_taxa()
# 
# colnames(data_fishbase) <- 
#   tolower(colnames(data_fishbase))
# data_fishbase$species <- gsub(" ", "_", data_fishbase$species)


## ---- load functions ----
melt.dist <- function(distmat, metric) {
  if(class(distmat)[1] == "dist") {distmat <- as.matrix(dist)}
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



read.genlight <- function(filters = "missind_callrate0.70_maf0.05",
                          level = "site",
                          removeless2ind = FALSE,
                          site2drop = NULL,
                          site2keep = NULL,
                          station2drop = NULL,
                          station2keep = NULL){
  
  # read genlight  
  genlight <- readRDS(paste0("intermediate/01_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, ".RDS"))
  
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
  if (level == "station"){
    genlight@pop <- 
      genlight@other[["ind.metrics"]][["site"]] %>%
      droplevels()
  }

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
  
  # remove populations with less than two individuals
  if (removeless2ind == TRUE){
    print(paste("removing", level, ":", names(which(table(genlight@pop) < 2))))
    genlight <- gl.drop.pop(genlight, pop.list = names(which(table(genlight@pop) < 2)), recalc = T, mono.rm = T)
  }

  # return
  return(genlight)  
}





# # check species database
# temp1 <-
#   data_species %>%
#   dplyr::filter(family == "Lutjanidae")
# 
# temp2 <-
#   data_species2 %>%
#   dplyr::filter(family == "Lutjanidae")
# 
# species1 <- temp1$species
# species2 <- temp2$species
# 
# species1[!(species1 %in% species2)]
# species2[!(species2 %in% species1)]



## ---- create arborescence ----
dir.create("intermediate/", showWarnings = F)
dir.create("intermediate/00_sampling_sites/", showWarnings = F)
dir.create("intermediate/01_genetic_diversity/", showWarnings = F)
dir.create("intermediate/02_species_diversity/", showWarnings = F)
dir.create("intermediate/03_distance_decay/", showWarnings = F)
dir.create("intermediate/05_re_Lesturgie/", showWarnings = F)

dir.create("results/", showWarnings = F)
dir.create("results/00_sampling_sites/", showWarnings = F)
dir.create("results/01_genetic_diversity/", showWarnings = F)
dir.create("results/02_species_diversity/", showWarnings = F)
dir.create("results/03_distance_decay/", showWarnings = F)
dir.create("results/04_continuity/", showWarnings = F)
dir.create("results/05_re_Lesturgie/", showWarnings = F)



## ---- load scripts ----
# 1. ...
# source("scripts/01_genetic_diversity/01_snp_fitering.R")




