# Wrapper script for the Etelis coruscans analyses

## ---- libraries ----
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
library(viridis)

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



## ---- data ----
# genetic data
data_samples <- read.csv("data/metadata_samples.csv")
data_sites <- read.csv("data/metadata_sites.csv")

data_samples <- 
  data_samples %>% 
  left_join(data_sites, by = c("station", "site")) %>%
  arrange(order)

# data_samples[,1:19] %>%
#   write.csv("intermediate/0_sampling_design/metadata_samples_ordered.csv", row.names = F, quote = F, na = "")

# species data
data_Etelis <- readRDS("data/Presence_data_Fishbase_Etelis_coruscans.RDS")

# taxonomy data
data_species <- read.csv("data/data_species_depth_range_teleo.csv")
# data_species2 <- read.csv("data/data_species.csv")
# data_fishtree <- read.csv("data/PFC_taxonomy.csv")
# data_fishbase <- rfishbase::load_taxa()
# 
# colnames(data_fishbase) <- 
#   tolower(colnames(data_fishbase))
# data_fishbase$species <- gsub(" ", "_", data_fishbase$species)


## ---- functions ----
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




## ---- arborescence ----
dir.create("intermediate/", showWarnings = F)
dir.create("intermediate/0_sampling_design/", showWarnings = F)
dir.create("intermediate/1_genetic_diversity/", showWarnings = F)
dir.create("intermediate/2_species_diversity/", showWarnings = F)
dir.create("intermediate/3_distance_decay/", showWarnings = F)
dir.create("intermediate/5_re_Lesturgie/", showWarnings = F)

dir.create("results/", showWarnings = F)
dir.create("results/0_sampling_design/", showWarnings = F)
dir.create("results/1_genetic_diversity/", showWarnings = F)
dir.create("results/2_species_diversity/", showWarnings = F)
dir.create("results/3_distance_decay/", showWarnings = F)
dir.create("results/4_continuity/", showWarnings = F)
dir.create("results/5_re_Lesturgie/", showWarnings = F)



## ---- scripts ----
# 1. ...
# source("scripts/1_genetic_diversity/11_snp_fitering.R")




