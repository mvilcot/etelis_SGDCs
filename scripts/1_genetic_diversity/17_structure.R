# library(fields)
# library(RColorBrewer)
# library(mapplots)



## ---- read SNPs dataset ----

# parameters
filters = "missind1_callrate0.70_maf0.05"
level = "station"
sites = "Hawaii"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = "Hawaii",
                station2drop = NULL,
                station2keep = NULL)


## ---- convert to LEA data format ----

# SNP presence/absence lfmm (package LEA, SilicoDArT)
gl2geno(genlight, 
        outfile = paste0("GenoLEA_Etelis_coruscans_ordered_", sites),
        outpath='./intermediate/1_genetic_diversity')


## ---- snmf ----
library(LEA)

# run snmf
obj.snmf <-
  LEA::snmf(paste0("intermediate/1_genetic_diversity/GenoLEA_Etelis_coruscans_ordered_", sites, ".geno"), 
       K = 1:10, alpha = 100, project = "new", repetitions = 1, entropy = TRUE)

# plot 
pdf(paste0("results/1_genetic_diversity/snmf_", sites, ".pdf"),
    height = 8, width = 11)
plot(obj.snmf, cex = 1.2, col = "lightblue", pch = 19)
for (K in 2:10){
  # get qmatrix
  qmatrix <- Q(obj.snmf, K = K, run = best)
  # get the cross-entropy of the 10 runs for K = 4
  ce = cross.entropy(obj.snmf, K = bestK)
  # select the run with the lowest cross-entropy for K = 4
  best = which.min(ce)
  # barplot
  barplot(t(qmatrix), col=RColorBrewer::brewer.pal(9,"Paired"), 
          border=NA, space=0, xlab="Individuals", 
          ylab="Admixture coefficients")
  # add legend
  for (i in 1:length(levels(genlight@pop))){
    axis(1, at=median(which(genlight@pop==levels(genlight@pop)[i])), labels=levels(genlight@pop)[i])}
}
dev.off()



## ---- tess3r ----
library(tess3r)



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
