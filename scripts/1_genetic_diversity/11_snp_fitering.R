
# ---- load ----
data_samplesRAW <- read_csv("data/metadata_samples_full.csv")


# ---- merge datasets ----
DART_Seamounts <- 
  read.csv("data/Report_DEtel22-6705_SNP_2_ordered_Seamounts.csv", header = F)

DART_Bowen <- 
  read.csv("data/Report_DEtel22-6705_SNP_2_ordered_Bowen.csv", header = F)

# merge
DART <- 
  DART_Bowen %>% 
  cbind(DART_Seamounts[, -c(1:18)])

# get individuals order from metadata
ind_order <- 
  match(data_samplesRAW$id, DART[7, ]) # DART 7th row stores individuals names

# reorder Dart
DART <- 
  DART[, c(1:18, ind_order)]

cat("nb individuals =", length(ind_order))

# save reordered dart SNP file
DART %>% 
  write.table("intermediate/1_genetic_diversity/Report_DEtel22-6705_SNP_2_ordered_full.csv",
              sep = ",", row.names = F, col.names = F, quote = F)


# ---- setup genlight ----

# read sorted DART file as genlight
gl <- 
  gl.read.dart(filename = "intermediate/1_genetic_diversity/Report_DEtel22-6705_SNP_2_ordered_full.csv",
               ind.metafile = "data/metadata_samples_full.csv") # takes the order of the DART columns for individuals

# add site as pop info
gl@pop <- 
  gl@other$ind.metrics[["site"]]

# save as RDS 
gl <- readRDS(file = "intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered.RDS")

# number of SNPs and individuals: 105144 SNPs, 369 individuals
gl



# ---- Graphic visualisation ----
gl <- 
  readRDS(file = "intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered2.RDS")

gl.report.bases(gl)
gl.report.callrate(gl, method = 'loc')
gl.report.callrate(gl, method = 'ind') # missing data by individuals
gl.report.maf(gl) # mininum allele frequency
gl.report.reproducibility(gl)
gl.report.heterozygosity(gl)
gl.report.rdepth(gl)
gl.report.secondaries(gl)

# get callrate by individual
options(max.print=2000)
t1 <- capture.output(gl.report.callrate(gl, method = "ind"))
t2 <- t1[38:length(t1)-1]
t3 <- as.data.frame(t2)
t4 <- 
  read.table(text=sub("^(\\S+)\\s+.*\\s+(\\S+)$", "\\1 \\2", t3$t2),
                 header=FALSE, stringsAsFactors= FALSE) %>% 
  as_tibble()
colnames(t4) <- t4[1,]
t4 <- t4[-1,]
t5 <- 
  t4 %>%
  dplyr::rename(id = "ind_name") %>% 
  dplyr::rename(site = "pop") %>% 
  dplyr::rename(callrate = "missing_data")

options(max.print=1000)

t5bis %>% 
  write_csv("results/1_genetic_diversity/gl_report_callrate_ind_Etelis_coruscans.csv")


# ---- Filtering ----

## 1 - Individual callrate ----
# ALSO REMOVE POPULATIONS IN WHICH LESS THAN TWO INDIVIDUALS
# Filter individuals called below 50%, identified with vcftools
# 84913 SNPs, 364 individuals
# Delete monomorphic loci here, that's why we have a different number of SNPs from vcftools pipeline
ind_to_remove <-
  t5 %>% 
  dplyr::filter(callrate < 0.5) %>% 
  dplyr::pull(id)
  
# remove c("ECO0934", "ECO0910", "ECO0912", "ECO0534", "ECO0535")
gl1 <- 
  gl %>% 
  gl.drop.ind(ind_to_remove, recalc = T, mono.rm = T) 
gl1

## 2 - Populations with only 1 indiviual left ----
# 84875 SNPs, 363 individuals
# "ECO1111" removed
gl2 <- 
  gl1 %>% 
  gl.drop.pop(pop.list = names(which(table(gl1@pop) < 2)), recalc = T, mono.rm = T)
gl2


## 3 - Callrate by site - 69304 SNPs ----
# Filter genotypes called below 70%
gl3 <- 
  gl2 %>% 
  gl.filter.callrate(threshold = 0.70)
gl3

## 4 - MAF - 21948 SNPs ----
# Filter genotypes with a allele frequency below 5%
gl4 <- 
  gl3 %>% 
  gl.filter.maf(threshold = 0.05)
gl4


# ---- relevel location factors ----

# relevel metadata
for (level in c("site", "station")){
  gl4$other$ind.metrics[[level]] <- 
    gl4$other$ind.metrics[[level]] %>% 
    ordered(levels = unique(data_stations[order(data_stations$order),][[level]])) %>% 
    droplevels()
}

# add site as pop info
gl4@pop <- 
  gl4@other$ind.metrics[["site"]]


# ---- subset metadata ----
# data_samples <-
#   data_samplesRAW %>%
#   dplyr::filter(id %in% gl4$ind.names)
# data_samples %>% write.csv("data/metadata_samples.csv",
#                               row.names = F, quote = T)



# ---- Save ----
filters <- "missind1_callrate0.70_maf0.05"

gl4 %>% 
  saveRDS(paste0("intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, ".RDS"))


# ---- Smearplot ----
gl4 %>%
  gl.smearplot(plot_colors = c("blue", "red", "green", "white"))

