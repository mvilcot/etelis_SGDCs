
library("dartR")
library("ggplot2")
library("wesanderson")


# ----------- Read dataset -----------------------------------------------------

## Set parameters
# filters <- "missind_callrate0.70_1stSNP_reprod1_rdepth_maf0.01"
filters <- "missind_callrate0.70_maf0.05"
filters <- "missind_callrate0.70_maf0.05_pcadapt"
sites <- "allsites"
sites <- "noSEYCH"
sites <- "Hawaii"

genlight <- readRDS(paste0("Intermediate/Genlight_DartSeq_Etelis_coruscans_grouped_", filters, "_", sites, ".RDS"))
genlight



# ----------- DAPC with sampling sites as prior --------------------------------

## DAPC - populations as prior
set.seed(999) # Setting a seed for a consistent result
dapc <- dapc(genlight)
saveRDS(dapc, paste0("Intermediate/DAPC_result_Etelis_coruscans_grouped_", filters, "_", sites, "_popprior.RDS"))

## Results
dapc <- readRDS(paste0("Intermediate/DAPC_result_Etelis_coruscans_grouped_", filters, "_", sites, "_popprior.RDS"))
print.dapc(dapc)
summary.dapc(dapc)
predict.dapc(dapc)

## Plots
png(paste0("Results/04_DAPC/DAPC_scatter_Etelis_coruscans_grouped_", filters, "_", sites, "_popprior.png"),
    height = 12, width = 12,
    units = 'in', res = 300)
scatter(dapc, scree.pca = F, legend = T, clabel = F, cleg = 0.8)
dev.off()

png(paste0("Results/04_DAPC/DAPC_compoplot_Etelis_coruscans_grouped_", filters, "_", sites, "_popprior.png"),
    height = 8, width = 20,
    units = 'in', res = 100)
compoplot(dapc, legend = T, show.lab = F, cleg=0.6)
dev.off()


## Plot ggplot - adapted from https://github.com/laurabenestan/TD_science_conservation

# Order labels of pop by longitude
dapc <- readRDS(paste0("Intermediate/DAPC_result_Etelis_coruscans_grouped_", filters, "_", sites, "_popprior.RDS"))
coord <- read.csv("Data/sites_coordinates.csv")
colnames(coord)[2] <- 'station'
metadata <- read.csv("Data/sample_metadata_grouped.csv")
pop <- merge(metadata, coord, by = "station", all.x = T)
colnames(pop)[2] <- 'INDIVIDUALS'

tab = as.data.frame(dapc$ind.coord)
tab$INDIVIDUALS <- row.names(tab)
dapc_geo <- merge(x = tab, y = pop, by=c("INDIVIDUALS"))

dapc_geo$pop <- factor(dapc_geo$station, levels=unique(dapc_geo[order(dapc_geo$order),]$station), ordered=TRUE)

# Plot
g <- ggplot(dapc_geo, aes(x=LD1, y=LD2, color=station)) +
  geom_point(size=3, alpha = 0.3) +
  scale_color_manual(values=wes_palette('Zissou1', nlevels(dapc_geo$pop), 'continuous')) +
  labs(x="DPC1") +
  labs(y="DPC2") +
  geom_text(
    label=dapc_geo$INDIVIDUALS,
    nudge_x = 0.25, nudge_y = 0.25,
    check_overlap = T, show.legend = FALSE
  )+
  # coord_fixed(ratio = 1) +
  theme_classic()
g
ggsave(paste0("Results/04_DAPC/DAPC_scatterperso_Etelis_coruscans_grouped_", filters, "_", sites, "_popprior.png"),
       height = 8, width = 14)



# ----------- TO DO!! - DAPC with clustering ---------------------------------
grp <- find.clusters(genlight, max.n.clust=20)
table(pop(genlight), grp$grp)
dapc1 <- dapc(genlight, grp$grp)
scatter(dapc1, scree.pca = F, legend = T, clabel = F, cleg = 0.8)
compoplot(dapc1, legend = T, show.lab = F, cleg=0.6)

saveRDS(dapc, paste0("Intermediate/DAPC_result_Etelis_coruscans_", filters, "_allsites_5groups.RDS"))




# ----------- TO DO!! - PCA --------------------------------------------------

PCA <- glPca(genlight, center = TRUE, scale = FALSE, nf = 10, loadings = TRUE, 
             alleleAsUnit = FALSE, useC = TRUE, parallel = require("parallel"),
             n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

DataPlot <- as.data.frame(PCA$scores)
DataPlot$pop <- genlight$pop
ggplot(DataPlot) +
  geom_point(aes(PC1, PC2, color=pop)) +
  theme_classic()





