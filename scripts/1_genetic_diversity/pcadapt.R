
# ---- PCADAPT ----
## ---- read Genlight into PCAdapt format ----

# parameters
sites = "allsites"
filters = "missind1_callrate0.70_maf0.05"
level = "site"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)




# transform to data frame
gdf <- as.data.frame(genlight)
genotype <- t(gdf)
gdf$pop <- genlight$pop # set population information

# transform to pcadapt format
library(pcadapt)
genotype_pca <- read.pcadapt(genotype)



## ---- PCAdapt ----

# first run
PCADAPT <- pcadapt(genotype_pca, K = 15)

# check number of PCs (=K) to retain
plotSCREE <- plot(PCADAPT, option = "screeplot") 
ggsave(paste0("results/1_genetic_diversity/PCAdapt_Screeplot_", filters, "_", sites, ".pdf"), 
       plotSCREE, height = 5, width = 8)

# rerun with good number of PCs
K <- 2 # K=1 or 2 with all sites
PCADAPT <- pcadapt(genotype_pca, K = K)


## ---- Graphic visualisation ----

# plot PC1-PC2
plotPCs <- plot(PCADAPT, option = "scores", i = 1, j = 2, pop = gdf$pop)
ggsave(paste0("results/1_genetic_diversity/PCAdapt_PC1PC2_", filters, "_", sites, "_K", K, ".pdf"), 
       plotPCs, height = 5, width = 8)

# Manhattan plot
plotMH <- plot(PCADAPT, option = "manhattan") 
ggsave(paste0("results/1_genetic_diversity/PCAdapt_Manhattan_", filters, "_", sites, "_K", K, ".pdf"), 
       plotMH, height = 5, width = 8)

# Check the expected uniform distribution of the p-values
plotQQ <- plot(PCADAPT, option = "qqplot", threshold = 0.1) 
ggsave(paste0("results/1_genetic_diversity/PCAdapt_QQplot_", filters, "_", sites, "_K", K, ".pdf"), 
       plotQQ, height = 5, width = 8)

# The excess of small p-values indicates the presence of outliers
plotDIST <- plot(PCADAPT, option = "stat.distribution") 
ggsave(paste0("results/1_genetic_diversity/PCAdapt_Distribution_", filters, "_", sites, "_K", K, ".pdf"), 
       plotDIST, height = 5, width = 8)

# presence of outliers is also visible when plotting a histogram of the test statistic Dj.
plotHIST <- ggplot() + xlab("p-values") +
  geom_histogram(aes(PCADAPT$pvalues), bins = 50, fill = "orange", colour = "black")
ggsave(paste0("results/1_genetic_diversity/PCAdapt_Histogram_", filters, "_", sites, "_K", K, ".pdf"), 
       plotHIST, height = 5, width = 8)



## ---- Choose a cutoff for outlier detection ----
library(qvalue)
# q-values
qval <- qvalue(PCADAPT$pvalues)$qvalues
alpha <- 0.1
outliersQ <- which(qval < alpha)
length(outliersQ) ## N = 3056 outliers

# Bonferroni correction (conservative)
padj <- p.adjust(PCADAPT$pvalues,method="bonferroni")
alpha <- 0.1
outliersB <- which(padj < alpha)
length(outliersB) ## N = 1744 outliers

# compare both methods
setdiff(outliersB, outliersQ) # in outliersB but not in outliersQ
setdiff(outliersQ, outliersB) # in outliersQ but not in outliersB

# keep outliers based on 10% Q-values
outliers <- outliersB

# extract corresponding SNP names and positions
outlier_loci <- data.frame(loci = rownames(genotype)[outliers], 
                           position = outliers)
write.table(outlier_loci, 
            file = paste0("intermediate/1_genetic_diversity/pcadapt_outliersB_loci_position_", sites, ".txt"), 
            sep = "\t", quote = F, row.names = F, col.names = F)

cat(length(outliers), "outlier loci")



## ---- Remove outlier loci from genlight ----
genlightADAPT <- gl.drop.loc(genlight, loc.list = outlier_loci$loci)
genlightADAPT 
# 18892 SNPs remaining

genlightADAPT %>% 
  saveRDS(paste0("intermediate/1_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, "_PCADAPTb_", sites, ".RDS"))





# ---- OUTFLANK ----
## ---- read Genlight into PCAdapt format ----

# parameters
sites = "allsites"
filters = "missind1_callrate0.70_maf0.05"
level = "site"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)

# run ----

temp <- gl.outflank(genlight, qthreshold = 0.01)

