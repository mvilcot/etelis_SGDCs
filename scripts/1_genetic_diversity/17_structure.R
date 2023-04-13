library(fields)
library(RColorBrewer)
library(mapplots)



## ---- read SNPs dataset ----

# parameters
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "noCocos"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)


## ---- convert to LEA data format ----

# SNP presence/absence lfmm (package LEA, SilicoDArT)
gl2geno(genlight, 
        outfile = paste0("Genlight_Etelis_coruscans_ordered_", filters),
        outpath='./intermediate/1_genetic_diversity')


## ---- run snmf ----

## With LEA
library(LEA)

# run snmf
obj.snmf <-
  LEA::snmf(paste0("intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, ".geno"), 
       K = 1:10, alpha = 100, project = "new", repetitions = 10, entropy = TRUE)

# get best K value
# plot cross-entropy criterion of all runs of the project
plot(obj.snmf, cex = 1.2, col = "lightblue", pch = 19)
# get the cross-entropy of the 10 runs for K = 4
ce = cross.entropy(obj.snmf, K = 2)
# select the run with the lowest cross-entropy for K = 4
best = which.min(ce)

# get Qmatrix
qmatrix <- Q(obj.snmf, K = 3)



## ---- plot ----

# plot qmatrix
par(mar=c(4,4,0.5,0.5))
barplot(t(qmatrix), col=RColorBrewer::brewer.pal(9,"Paired"), 
        border=NA, space=0, xlab="Individuals", 
        ylab="Admixture coefficients")
for (i in 1:length(levels(genlight@pop))){
  axis(1, at=median(which(genlight@pop==levels(genlight@pop)[i])), labels=levels(genlight@pop)[i])}





## ---- STRUCTURE ----
library(strataG)
structure <- gl.run.structure(
  genlight,
  k.range = 2:5,
  num.k.rep = 3,
  # exec = "C:/Users/User/Downloads/structure_windows_console/console/structure.exe",
  plot.out = TRUE,
  save2tmp = FALSE,
  verbose = NULL
)
