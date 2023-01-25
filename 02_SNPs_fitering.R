
library("dartR") # gl.read.dart
library("ggtern") # hwe
library("ggplot2")


# ----------- Read DART file into genlight -------------------------------------

## Read DART file and metadata
DART <- read.csv("Data/DartSeq/Report_DEtel22-6705_SNP_2.csv")
metadata <- read.csv("Data/sample_metadata_grouped.csv")
# metadata <- read.csv("Data/sample_metadata.csv")

## Reorder individuals according to populations
DART <- DART[,c(1:18, match(metadata$id, DART[6,]))]
## !!!!!!WARNING!!!!!! Don't rewrite this csv, lines 1 and 6 must be modified by hand
# write.csv(DART, "Intermediate/Report_DEtel22-6705_SNP_2_grouped.csv",
          # quote = F, row.names = F)

## Read sorted DART file as GENLIGHT
gl <- dartR::gl.read.dart(filename="Intermediate/Report_DEtel22-6705_SNP_2_grouped.csv",
                          ind.metafile="Data/sample_metadata_grouped.csv")

## Save as .RDS 
## I prefer this format compared to Rdata because you can assign to another variable name
saveRDS(gl, file = "Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped.RDS")

## Observe number of SNPs and individuals: 105144 SNPs, 369 individuals
gl




# ----------- Graphic visualisation ---------------------------------------------
gl <- readRDS(file = "Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped.RDS")

gl.report.bases(gl)
gl.report.callrate(gl, method = 'loc')
gl.report.callrate(gl, method = 'ind') # missing data by individuals
gl.report.maf(gl) # mininum allele frequency
gl.report.reproducibility(gl)
gl.report.heterozygosity(gl)
gl.report.rdepth(gl)




# ----------- Filtering ---------------------------------------------------------
## 1. Individual callrate
# Filter individuals called below 50%, identified with vcftools
# 84913 SNPs, 364 individuals
# Delete monomorphic loci here, that's why we hane a different number of SNPs from vcftools pipeline
gl1 <- gl.drop.ind(gl, c("ECO0934", "ECO0910", "ECO0912", "ECO0534", "ECO0535"), 
                   recalc = F, mono.rm = T) 
gl1

## 2. Callrate by site - 69404 SNPs
# Filter genotypes called below 70%
gl2 <- gl.filter.callrate(gl1, threshold = 0.70)
gl2

## 3 - MAF - 21948 SNPs
# Filter genotypes with a allele frequency below 5%
gl3 <- gl.filter.maf(gl2, threshold = 0.05)
gl3

## Save filtered SNPs data
# filters <- "missind_callrate0.70_1stSNP_reprod1_rdepth_maf0.01"
filters <- "missind_callrate0.70_maf0.05"
saveRDS(gl3, paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_allsites.RDS"))




# ----------- Generate subset without seychelles --------------------------------
gl3_noSEYCH <- gl.drop.pop(gl3, pop.list = c("Seychelles"), recalc = T, mono.rm = T)
saveRDS(gl3_noSEYCH, paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_noSEYCH.RDS"))




# ----------- Smearplot ---------------------------------------------------------
gl.smearplot(gl3, plot_colors = c("blue", "red", "green", "white"))




# ----------- OLD filters --------------------------------------------------------

# ## 2b - keep 1st SNP of each fragment - 44392 SNPs
# gl2b <- gl.filter.secondaries(gl2, method = 'best')
# 
# ## 2c - Reproducibility - 29826 SNPs
# gl2c <- gl.filter.reproducibility(gl2b, threshold = 1)
# 
# ## 2d - Maximum read depth - 26194 SNPs
# gl.report.rdepth(gl2c)
# d <- 16.33859 # mean read depth
# gl2d <- gl.filter.rdepth(gl2c, upper = d+4*sqrt(d), lower = 0)
# gl.report.rdepth(gl2d)
# 
# ## !!! departure from Hardy-Weinberg equilibrium
# gl.report.hwe(gl3, subset = 'all') # all samples together
# gl.report.hwe(gl3, subset = 'each') # each population sparatedly




