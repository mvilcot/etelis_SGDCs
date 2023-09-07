
# ---- read SNPs dataset ----

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



# ---- PCA ----

## ---- run ----
# PCA <- adegenet::glPca(genlight, center = TRUE, scale = FALSE, nf = 10, loadings = TRUE,
#                        alleleAsUnit = FALSE, useC = TRUE, parallel = require("parallel"),
#                        n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

# PCA <- dartR::gl.pcoa(genlight)
# saveRDS(PCA, paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, ".RDS"))
PCA <- readRDS(paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, ".RDS"))



## ---- plot auto ----
png(paste0("results/1_genetic_diversity/PCA_auto_", filters, "_", sites, "_", level, ".png"),
    height = 6, width = 7, units = "in", res = 300)
dartR::gl.pcoa.plot(PCA, genlight)
dev.off()


## ---- plot perso ----
# get dapc values and metadata
pca_scores <- 
  as.data.frame(PCA$scores) %>% 
  rownames_to_column("id") %>% 
  dplyr::left_join(data_samples, by = "id")

# setup axis
percent_explained <- PCA$eig / sum(PCA$eig) * 100
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PCA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCA Axis 2 ({pretty_pe[2]}%)"))

# save  
gg_pca <- 
  ggplot(pca_scores, aes(x=PC1, y=PC2, color=.data[[level]])) +
  geom_point(size = 2, alpha = 0.6) +
  labs(x=labels[1], y=labels[2]) +
  # scale_color_manual(values = color_perso) +
  theme_classic()
ggsave(paste0("results/1_genetic_diversity/PCA_perso_", filters, "_", sites, "_", level, ".png"),
       gg_pca,
       height = 6, width = 8)


