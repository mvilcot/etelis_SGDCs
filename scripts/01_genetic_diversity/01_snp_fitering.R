

## ---- read DART file into genlight ----

# read DART output file
DART <- read.csv("data/Report_DEtel22-6705_SNP_2.csv", header = F)

# reorder individuals (columns) according to metadata
ind_order <- match(data_samples$id, DART[7, ])
print(length(ind_order))
DART <- DART[, c(1:18, ind_order)]

# save reordered dart SNP file
DART %>% 
  write.table("intermediate/01_genetic_diversity/Report_DEtel22-6705_SNP_2_ordered.csv",
              sep = ",", row.names = F, col.names = F)

# read sorted DART file as genlight
gl <- gl.read.dart(filename="intermediate/01_genetic_diversity/Report_DEtel22-6705_SNP_2_ordered.csv",
                   ind.metafile="data/metadata_samples.csv")

# add site as pop info
gl@pop <- gl@other[["ind.metrics"]][["site"]]


# save as RDS 
gl %>% 
  saveRDS(file = "intermediate/01_genetic_diversity/Genlight_Etelis_coruscans_ordered.RDS")

# number of SNPs and individuals: 105144 SNPs, 369 individuals
print(gl)




## ---- Graphic visualisation ----
gl <- readRDS(file = "intermediate/01_genetic_diversity/Genlight_Etelis_coruscans_ordered.RDS")

gl.report.bases(gl)
gl.report.callrate(gl, method = 'loc')
gl.report.callrate(gl, method = 'ind') # missing data by individuals
gl.report.maf(gl) # mininum allele frequency
gl.report.reproducibility(gl)
gl.report.heterozygosity(gl)
gl.report.rdepth(gl)




## ---- Filtering ----
# 1. Individual callrate
# Filter individuals called below 50%, identified with vcftools
# 84913 SNPs, 364 individuals
# Delete monomorphic loci here, that's why we hane a different number of SNPs from vcftools pipeline
gl1 <- gl.drop.ind(gl, c("ECO0934", "ECO0910", "ECO0912", "ECO0534", "ECO0535"), 
                   recalc = F, mono.rm = T) 
gl1

# 2. Callrate by site - 69404 SNPs
# Filter genotypes called below 70%
gl2 <- gl.filter.callrate(gl1, threshold = 0.70)
gl2

# 3 - MAF - 21948 SNPs
# Filter genotypes with a allele frequency below 5%
gl3 <- gl.filter.maf(gl2, threshold = 0.05)
gl3

# Save filtered SNPs data
filters <- "missind_callrate0.70_maf0.05"
gl3 %>% 
  saveRDS(paste0("intermediate/01_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, ".RDS"))



## ---- Smearplot ----
gl.smearplot(gl3, plot_colors = c("blue", "red", "green", "white"))



