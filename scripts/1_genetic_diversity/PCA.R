
# read SNPs dataset ----

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

genlight


# run PCA ----
## perform PCA
PCA <- glPca(genlight, center = TRUE, scale = FALSE, nf = 10, loadings = TRUE, 
             alleleAsUnit = FALSE, useC = TRUE, parallel = FALSE,
             n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

DataPlot <- as.data.frame(PCA$scores)
DataPlot$pop <- genlight$pop


## assign rad IDs
pca_xy <- 
  pca_run$li %>% 
  rownames_to_column("rad_id")

pca_xy <- 
  data_meta %>% 
  dplyr::right_join(pca_xy, by = "rad_id")

## relevel factor
pca_xy$sample_site <- factor(pca_xy$sample_site, levels = names(site_colours))

## plot
gg <-
  ggplot(pca_xy, aes(x=Axis1, y=Axis2, color = sample_site, label = rad_id)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = site_colours,
                     name = "Sampling site",
                     breaks = names(site_colours)) +
  labs(color = "Sample site") +
  theme_classic() +
  ggtitle(sp)

# ggplotly(gg)

ggsave(paste0("results/PCA/", ocean, "_PCA_allsamples_allSNPs_", sp, "_samples_excluded.png"),
       plot = gg,
       height = 8, width = 8)


