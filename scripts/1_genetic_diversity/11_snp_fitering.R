
## ---- read Dart file into genlight ----
# read DART output file
DART <- 
  read.csv("data/Report_DEtel22-6705_SNP_2.csv", header = F)


## ---- reorder individuals ----
# reorder Dart columns (individuals) according to metadata order 


# get individuals order from metadata
ind_order <- 
  match(data_samples$id, DART[7, ]) # DART 7th row stores individuals names

# reorder Dart
DART <- 
  DART[, c(1:18, ind_order)]

cat("nb individuals =", length(ind_order))

# save reordered dart SNP file
DART %>% 
  write.table("intermediate/1_genetic_diversity/Report_DEtel22-6705_SNP_2_ordered.csv",
              sep = ",", row.names = F, col.names = F)


## ---- setup genlight ----

# read sorted DART file as genlight
gl <- 
  gl.read.dart(filename = "intermediate/1_genetic_diversity/Report_DEtel22-6705_SNP_2_ordered.csv",
               ind.metafile = "data/metadata_samples.csv") # takes the order of the DART columns for individuals

# add site as pop info
gl@pop <- 
  gl@other$ind.metrics[["site"]]

# save as RDS 
gl %>% 
  saveRDS(file = "intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered.RDS")

# number of SNPs and individuals: 105144 SNPs, 369 individuals
gl



## ---- Graphic visualisation ----
gl <- 
  readRDS(file = "intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered.RDS")

# gl.report.bases(gl)
# gl.report.callrate(gl, method = 'loc')
# gl.report.callrate(gl, method = 'ind') # missing data by individuals
# gl.report.maf(gl) # mininum allele frequency
# gl.report.reproducibility(gl)
# gl.report.heterozygosity(gl)
# gl.report.rdepth(gl)

# get callrate by individual
options(max.print=2000)
t1 <- capture.output(gl.report.callrate(gl, method = "ind"))
t2 <- t1[38:length(t1)-1]
t3 <- as.data.frame(t2)
t4 <- read.table(text=sub("^(\\S+)\\s+.*\\s+(\\S+)$", "\\1 \\2", t3$t2),
                 header=FALSE, stringsAsFactors= FALSE)
colnames(t4) <- t4[1,]
t4 <- t4[-1,]
t4 %>%
  rename(id = ind_name)

options(max.print=1000)

write.csv(t4, "results/1_genetic_diversity/gl_report_callrate_ind_Etelis_coruscans.csv",
          quote = F, row.names = F)


## ---- Filtering ----

# 1 - Individual callrate
# !!!! ALSO REMOVE POPULATIONS IN WHICH LESS THAN TWO INDIVIDUALS
# Filter individuals called below 50%, identified with vcftools
# 84913 SNPs, 364 individuals
# Delete monomorphic loci here, that's why we have a different number of SNPs from vcftools pipeline
gl1 <- 
  gl %>% 
  gl.drop.ind(c("ECO0934", "ECO0910", "ECO0912", "ECO0534", "ECO0535"), recalc = T, mono.rm = T) 
gl1

# 2 - Sites with only 1 indiviual left
# 84875 SNPs, 363 individuals
gl2 <- 
  gl1 %>% 
  gl.drop.pop(pop.list = names(which(table(gl1@pop) < 2)), recalc = T, mono.rm = T)
gl2


# 3 - Callrate by site - 69304 SNPs
# Filter genotypes called below 70%
gl3 <- 
  gl2 %>% 
  gl.filter.callrate(threshold = 0.70)
gl3

# 4 - MAF - 21948 SNPs
# Filter genotypes with a allele frequency below 5%
gl4 <- 
  gl3 %>% 
  gl.filter.maf(threshold = 0.05)
gl4


## ---- relevel location factors ----

# relevel metadata
for (level in c("site", "station")){
  gl4$other$ind.metrics[[level]] <- 
    gl4$other$ind.metrics[[level]] %>% 
    ordered(levels = unique(data_sites[order(data_sites$order),][[level]])) %>% 
    droplevels()
}

# add site as pop info
gl4@pop <- 
  gl4@other$ind.metrics[["site"]]


# ## ---- subset metadata ----
# data_samples <-
#   data_samples %>%
#   dplyr::filter(id %in% gl4$ind.names)
# data_samples %>% write.csv("intermediate/0_sampling_design/metadata_samples_subset.csv",
#                               row.names = F, quote = T)
# 
# data_sites <-
#   data_sites %>%
#   dplyr::filter(station %in% data_samples$station)
# data_sites %>% write.csv("intermediate/0_sampling_design/metadata_sites_subset.csv",
#                          row.names = F, quote = T)



## ---- Save ----
filters <- "missind1_callrate0.70_maf0.05"

gl4 %>% 
  writeRDS(paste0("intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, ".RDS"))


## ---- Smearplot ----
gl4 %>% 
  gl.smearplot(plot_colors = c("blue", "red", "green", "white"))

