# Wrapper script for the Etelis coruscans analyses

# ---- libraries ----
# scripts
library(tidyverse)      # the beautiful beautiful tidyverse
library(pracma)         # findintervals
library(glue)           # use braces for PXA axis
# library(reshape)        # melt
library(ecodist)        # MRM
library(dotwhisker)     # dotwhisket plot

# plot
library(ggplot2)        # plots
# library(wesanderson)    # palette
library(ggh4x)
library(patchwork)      # easy multiple plot
library(gridExtra)      # easy multiple plot
# library(paletteer)      # palette
library(viridis)        # palette

# spatial
library(terra)          # raster-type package
library(sdmpredictors)  # bio-oracle portal
library(geodist)        # geographic distance
library(marmap)         # bathymetry data
library(gdistance)
library(sf)
library(mapview)

# taxonomy
library(fishtree)       # fish tree of life 
library(betapart)       # jaccard diversity
library(ape)

# genetics
library(adegenet)
library(dartR)          # SNPs filtering
library(pcadapt)
# library(qvalue)
# library(OutFLANK)
library(hierfstat)
library(mmod)
library(pegas)          # genind2loci
library(jacpop)         # generate_pw_jaccard

# phylogenetics
library(picante)        # pd, mpd
# library(ggtree)         # beautiful tree plot
library(phytools)       # midpoint.root
library(ape)



# ---- data ----

## genetic data
# data_samples <- read.csv("data/metadata_samples.csv")
# data_sites <- read.csv("data/metadata_sites.csv")
# 
# data_samples <-
#   data_samples %>%
#   left_join(data_sites, by = c("station", "site")) %>%
#   arrange(order)

data_samples <- read.csv("intermediate/0_sampling_design/metadata_samples_subset.csv")

for (level in c("site", "station")){ # order levels
  data_samples[[level]] <- 
    data_samples[[level]] %>% 
    ordered(levels = unique(data_samples[order(data_samples$order),][[level]])) %>% 
    droplevels()
}

# data_samples[,1:19] %>%
#   write.csv("intermediate/0_sampling_design/metadata_samples_ordered.csv", row.names = F, quote = F, na = "")



## ---- presence data ----
data_Etelis <- readRDS("data/Presence_data_Fishbase_Etelis_coruscans.RDS")



## ---- taxonomy data ----
data_fishbase <- rfishbase::load_taxa()
# data_species2 <- read.csv("data/data_species.csv")
# data_fishtree <- read.csv("data/PFC_taxonomy.csv")

data_fishbase$Species <- gsub(" ", "_", data_fishbase$Species)
# colnames(data_fishbase) <- tolower(colnames(data_fishbase))


# species_list <- 
#   data_fishbase %>% 
#   dplyr::filter(Family == "Lutjanidae") %>% 
#   pull(Species)
# 
# temp <- 
#   rfishbase::ecology(species_list)



## ---- depth data ----
data_depth <- read_csv("data/data_species_depth_range_teleo.csv")



## ---- spatial data ----
data_sites <- read_csv("intermediate/0_sampling_design/metadata_sites_subset.csv")

# relevel
for (level in c("site", "station")){
  data_sites[[level]] <- 
    data_sites[[level]] %>% 
    ordered(levels = unique(data_sites[order(data_sites$order),][[level]])) %>% 
    droplevels()
}

# get mean coordinates by location
coord_site <- 
  data_sites %>% 
  group_by(site) %>% 
  summarise(longitude=mean(Longitude_approx),
            latitude=mean(Latitude_approx)) %>% 
  arrange(factor(site, levels(data_samples$site))) %>%
  column_to_rownames("site")




# ---- personalised plot ----
level = "site"
color_perso <- c(viridis(length(levels(data_samples[[level]]))))
names(color_perso) <- levels(data_samples[[level]])



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
shift.lonst = function(x) {
  geom = st_geometry(x)
  st_geometry(x) = st_sfc(
    lapply(seq_along(geom), function(i) {
      geom[[i]][1] = ifelse(geom[[i]][1] < 0, geom[[i]][1] + 360, geom[[i]][1])
      return(geom[[i]])
    })
    , crs = st_crs(geom)
  )

  longitudecolumn <- grep("Longitude|longitude", colnames(x))
  x[[longitudecolumn]] <- sapply(x[[longitudecolumn]], function(i){ifelse(i < 0, i + 360, i)})
  return(x)
}

shift.lon <- function(df) {
  longitudecolumn <- grep("Longitude|longitude", colnames(df))
  df[[longitudecolumn]] <- sapply(df[[longitudecolumn]], function(i){ifelse(i < 0, i + 360, i)})
  return(df)
}



## remove location from distance matrix
loc = "Seychelles"

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

mat.subset(dist_mat[[1]], "Seychelles")







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




