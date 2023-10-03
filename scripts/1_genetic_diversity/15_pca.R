

# ---- PCA all sites ----

## load data ----

# parameters
sites = "noCocos"
filters = "missind1_callrate0.70_maf0.05"
level = "site"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)

genlight



## run PCA ----
# PCA <- adegenet::glPca(genlight, center = TRUE, scale = FALSE, nf = 10, loadings = TRUE,
#                        alleleAsUnit = FALSE, useC = TRUE, parallel = require("parallel"),
#                        n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

# PCA <- dartR::gl.pcoa(genlight)
# saveRDS(PCA, paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, ".RDS"))
PCA <- readRDS(paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, ".RDS"))



## plot auto ----
# png(paste0("results/1_genetic_diversity/PCA_auto_", filters, "_", sites, "_", level, ".png"),
#     height = 6, width = 7, units = "in", res = 300)
# dartR::gl.pcoa.plot(PCA, genlight)
# dev.off()


## plot perso ----
# get dapc values and metadata
pca_scores <- 
  as.data.frame(PCA$scores) %>% 
  rownames_to_column("id") %>% 
  dplyr::left_join(data_samples, by = "id")

# setup axis
percent_explained <- PCA$eig / sum(PCA$eig) * 100
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PC1 ({pretty_pe[1]}%)"),
            glue("PC2 ({pretty_pe[2]}%)"))

# save  
gg1 <- 
  ggplot(pca_scores, aes(x=PC1, y=PC2, color=.data[[level]], label = id)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(x=labels[1], y=labels[2]) +
  scale_color_manual(values = color_perso) +
  theme_light() +
  labs(tag = "(a)") 

plotly::ggplotly(gg1)


# ---- PCA no Seychelles ----

## load data ----

# parameters
sites = "noCocos_noSeychelles"
filters = "missind1_callrate0.70_maf0.05"
level = "site"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = "Seychelles",
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)

genlight



## run PCA ----
# PCA <- adegenet::glPca(genlight, center = TRUE, scale = FALSE, nf = 10, loadings = TRUE,
#                        alleleAsUnit = FALSE, useC = TRUE, parallel = require("parallel"),
#                        n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

# PCA <- dartR::gl.pcoa(genlight)
# saveRDS(PCA, paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, ".RDS"))
PCA <- readRDS(paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, ".RDS"))



## plot auto ----
# png(paste0("results/1_genetic_diversity/PCA_auto_", filters, "_", sites, "_", level, ".png"),
#     height = 6, width = 7, units = "in", res = 300)
# dartR::gl.pcoa.plot(PCA, genlight)
# dev.off()


## plot perso ----
# get dapc values and metadata
pca_scores <- 
  as.data.frame(PCA$scores) %>% 
  rownames_to_column("id") %>% 
  dplyr::left_join(data_samples, by = "id")

# setup axis
percent_explained <- PCA$eig / sum(PCA$eig) * 100
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PC1 ({pretty_pe[1]}%)"),
            glue("PC2 ({pretty_pe[2]}%)"))

# save  
gg2 <- 
  ggplot(pca_scores, aes(x=PC1, y=PC2, color=.data[[level]], label = id)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(x=labels[1], y=labels[2]) +
  scale_color_manual(values = color_perso) +
  theme_light() +
  theme(legend.position = "none") +
  labs(tag = "(b)") 

plotly::ggplotly(gg2)





# ---- Save plots together ----

gg1 + gg2 + plot_layout(guides = "collect") 
ggsave(paste0("results/1_genetic_diversity/PCA_perso_", filters, "_", level, "_2.png"),
       height = 5, width = 12)



# ---- PCA Hawaii ----

## load data ----

# parameters
sites = "Hawaii"
filters = "missind1_callrate0.70_maf0.05"
level = "station"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = "Hawaii",
                station2drop = NULL,
                station2keep = NULL)

genlight



## run PCA ----
# PCA <- adegenet::glPca(genlight, center = TRUE, scale = FALSE, nf = 10, loadings = TRUE,
#                        alleleAsUnit = FALSE, useC = TRUE, parallel = require("parallel"),
#                        n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

# PCA <- dartR::gl.pcoa(genlight)
# saveRDS(PCA, paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, ".RDS"))
PCA <- readRDS(paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, ".RDS"))



## plot auto ----
# png(paste0("results/1_genetic_diversity/PCA_auto_", filters, "_", sites, "_", level, ".png"),
#     height = 6, width = 7, units = "in", res = 300)
# dartR::gl.pcoa.plot(PCA, genlight)
# dev.off()


## plot perso ----
# get dapc values and metadata
pca_scores <- 
  as.data.frame(PCA$scores) %>% 
  rownames_to_column("id") %>% 
  dplyr::left_join(data_samples, by = "id")

# setup axis
percent_explained <- PCA$eig / sum(PCA$eig) * 100
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PC1 ({pretty_pe[1]}%)"),
            glue("PC2 ({pretty_pe[2]}%)"))

# save  
gg <- 
  ggplot(pca_scores, aes(x=PC1, y=PC2, color=.data[[level]], label = id)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(x=labels[1], y=labels[2]) +
  theme_light()

plotly::ggplotly(gg)

ggsave(paste0("results/1_genetic_diversity/PCA_perso_", filters, "_", level, "_", sites, ".png"),
       height = 5, width = 7)

