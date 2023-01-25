
library("pcadapt")
library("qvalue")
library("OutFLANK")
library("dartR")


# ----------- Read Genlight into PCAdapt format --------------------------------
## Set parameters
filters <- "missind_callrate0.70_maf0.05"
sites <- "allsites" 

## Read genlight
genlight <- readRDS(paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_", sites, ".RDS"))
genlight

## Transform Genlight as data frame
df <- as.data.frame(genlight)
genotype <- t(df)

## Set population information
df$pop <- genlight$pop

## Get into pcadapt format
genotype_pca <- read.pcadapt(genotype)



# ----------- Run PCAdapt ------------------------------------------------------

PCADAPT <- pcadapt(genotype_pca, K = 15)

## Check number of PCs (=K) to retain
plotSCREE <- plot(PCADAPT, option = "screeplot") 
ggsave(paste0("Results/03_PCAdapt/PCAdapt_", filters, "_", sites, "_Screeplot.pdf"), 
       plotSCREE, height = 5, width = 8)


## Rerun with good number of PCs
K <- 2 # K=1 or 2 with all sites, K=6 without Seychelles
PCADAPT <- pcadapt(genotype_pca, K = K)


# ----------- Graphic visualisation --------------------------------------------

## plot PC1-PC2
plotPCs <- plot(PCADAPT, option = "scores", i = 1, j = 2, pop = df$pop)
ggsave(paste0("Results/03_PCAdapt/PCAdapt_", filters, "_", sites, "_K", K, "_PC1PC2.pdf"), 
       plotPCs, height = 5, width = 8)

## Manhattan plot
plotMH <- plot(PCADAPT, option = "manhattan") 
ggsave(paste0("Results/03_PCAdapt/PCAdapt_", filters, "_", sites, "_K", K, "_Manhattan.pdf"), 
       plotMH, height = 5, width = 8)

## Check the expected uniform distribution of the p-values
plotQQ <- plot(PCADAPT, option = "qqplot", threshold = 0.1) 
ggsave(paste0("Results/03_PCAdapt/PCAdapt_", filters, "_", sites, "_K", K, "_QQplot.pdf"), 
       plotQQ, height = 5, width = 8)

## The excess of small p-values indicates the presence of outliers
plotDIST <- plot(PCADAPT, option = "stat.distribution") 
ggsave(paste0("Results/03_PCAdapt/PCAdapt_", filters, "_", sites, "_K", K, "_Distribution.pdf"), 
       plotDIST, height = 5, width = 8)

## presence of outliers is also visible when plotting a histogram of the test statistic Dj.
plotHIST <- ggplot() + xlab("p-values") +
  geom_histogram(aes(PCADAPT$pvalues), bins = 50, fill = "orange", colour = "black")
ggsave(paste0("Results/03_PCAdapt/PCAdapt_", filters, "_", sites, "_K", K, "_Histogram.pdf"), 
       plotHIST, height = 5, width = 8)



# ----------- Choose a cutoff for outlier detection ----------------------------
## q-values
qval <- qvalue(PCADAPT$pvalues)$qvalues
alpha <- 0.1
outliersQ <- which(qval < alpha)
length(outliersQ) ## N = 3056 outliers

## Bonferroni correction (conservative)
padj <- p.adjust(PCADAPT$pvalues,method="bonferroni")
alpha <- 0.1
outliersB <- which(padj < alpha)
length(outliersB) ## N = 1744 outliers

## compare both methods
setdiff(outliersB, outliersQ) # in outliersB but not in outliersQ
setdiff(outliersQ, outliersB) # in outliersQ but not in outliersB

## keep outliers based on 10% Q-values
outliers <- outliersQ 

## extract corresponding SNP names and positions
outlier_loci <- data.frame(loci=rownames(genotype)[outliers], position=outliers)
write.table(outlier_loci, file = "Intermediate/pcadapt_outliers_loci_position.txt", sep = "\t", quote = F, row.names = F, col.names = F)
print(paste0(length(outliers), " outlier loci"), quote = 0)



# ----------- Remove outlier loci from genlight - ------------------------------
genlight <- gl.drop.loc(genlight, loc.list = outlier_loci$loci)
genlight # 18892 SNPs remaining
 
filters <- "missind_callrate0.70_maf0.05_pcadapt"
saveRDS(genlight, paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_", sites, ".RDS"))


# ----------- Generate subset without seychelles -------------------------------
genlight <- readRDS(paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_pcadapt_allsites.RDS"))
genlight_noSEYCH <- gl.drop.pop(genlight, pop.list = c("Seychelles"), recalc = T, mono.rm = T)
genlight_noSEYCH
saveRDS(genlight_noSEYCH, paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_noSEYCH.RDS"))

# ----------- Hawaii ----------------------------------------------------------
genlight_Hawaii <- gl.keep.pop(genlight, pop.list = c("Hawaii"), recalc = T, mono.rm = T)
genlight_Hawaii@pop <- genlight_Hawaii@other[["ind.metrics"]][["station"]]
genlight_Hawaii@pop <- droplevels(genlight_Hawaii@pop)
saveRDS(genlight_Hawaii, paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_Hawaii.RDS"))

