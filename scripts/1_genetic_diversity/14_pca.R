
# ---- PCA all sites ----

## load data ----

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

genlight


## run PCA ----
npca <- length(genlight$ind.names) - 1

PCA <- adegenet::glPca(genlight, nf = npca)
PCA %>% saveRDS(paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, "_", npca, "PCA_PCADAPTb.RDS"))
# PCA <- readRDS(paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, "_", npca, "PCA.RDS"))



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
  scale_color_manual('', values = color_perso, labels = LABELS) +
  theme_light()

plotly::ggplotly(gg1)

#### {FIGURE S2} ####
ggsave(paste0("results/1_genetic_diversity/_S2_PCA_perso_", filters, "_", level, "_", sites, ".png"),
       gg1,
       height = 5, width = 7, dpi = 500)




# ---- PCA no Seychelles ----

## load data ----

# parameters
sites = "noSeychelles"
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
npca <- length(genlight$ind.names) - 1

PCA <- adegenet::glPca(genlight, nf = npca)
PCA %>% saveRDS(paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, "_", npca, "PCA.RDS"))
# PCA <- readRDS(paste0("intermediate/1_genetic_diversity/PCA_output_", filters, "_", sites, "_", level, "_", npca, "PCA.RDS"))




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
  theme(legend.position = "none")

plotly::ggplotly(gg2)





# ---- Save plots together ----

gg1 + gg2 + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave(paste0("results/1_genetic_diversity/PCA_perso_", filters, "_", level, ".png"),
       height = 5, width = 12)


