# Wrapper script for the Etelis coruscans analyses

## ---- load libraries ----
# scripts
library(tidyverse)      # the beautiful beautiful tidyverse
library(pracma)         # findintervals

# plot
library(ggplot2)        # plots
library(wesanderson)    # palette
library(ggh4x)
library(paletteer)      # palette
library(patchwork)

# spatial
library(terra)          # raster-type package
library(sdmpredictors)  # bio-oracle portal

# taxonomy
library(fishtree)       # fish tree of life 
library(betapart)       # jaccard diversity

# genetics
# library(ggtern) # hwe
library(adegenet)
library(dartR)          # SNPs filtering
library(pcadapt)
library(qvalue)
library(OutFLANK)


## ---- load scripts ----
# load in all the sample data that is in common use
source("scripts/01_load.R")


## ---- genetic diversity ----
# 1. ...
# source("scripts/01_genetic_diversity/01_snp_fitering.R")
