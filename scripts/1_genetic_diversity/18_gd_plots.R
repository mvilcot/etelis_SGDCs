
# ---- all sites ----
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "allsites"

kmin <- 1
kmax <- 9
nrep <- 5
k <- 2
nloci <- 21921

## load ----
dapc <- readRDS(paste0("intermediate/1_genetic_diversity/DAPC_popprior_output_", filters, "_", sites, "_", level, "_53PCA_9DA.RDS"))
structure <- readRDS(paste0("intermediate/1_genetic_diversity/STRUCTURE_", nloci, "loci_", sites, "_K", kmin, "-", kmax, "_r", nrep, ".RDS"))

scatter(dapc, scree.pca=TRUE)


##  scatter plot ----
# Plot ggplot - adapted from https://github.com/laurabenestan/TD_science_conservation

# setup axis
percent_explained <- dapc$eig / sum(dapc$eig) * 100
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("LD1 ({pretty_pe[1]}%)"),
            glue("LD2 ({pretty_pe[2]}%)"))


# get dapc values and metadata
dapc_geo <- 
  as.data.frame(dapc$ind.coord) %>% 
  rownames_to_column("id") %>% 
  dplyr::left_join(data_samples, by = "id") %>% 
  cbind(grp = dapc$grp)

# reorder labels by longitude
dapc_geo$grp <- 
  factor(dapc_geo$grp, 
         levels = unique(dapc_geo[order(dapc_geo$order),][[level]]))


# Plot
gg_dapc1 <- ggplot(dapc_geo, aes(x=LD1, y=LD2, color=grp)) +
  geom_vline(xintercept = 0, color = "darkgrey") +
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_point(size = 3, alpha = 0.5) +
  stat_ellipse(level = 0.67, show.legend = FALSE) +
  scale_color_manual('Site', values = color_perso, labels = LABELS) +
  labs(x=labels[1], y=labels[2]) +
  theme_light()
gg_dapc1



## compolot ----

# create an object with membership probabilities
qmat <- gl.plot.structure(structure, K=k, colors_clusters = viridis(k))

### preparing data ----
qmat_plot_large <- 
  as.data.frame(qmat[[1]]) %>% 
  dplyr::select(c("orig.pop", "Label", contains("cluster"))) 
qmat_plot_long <- 
  qmat_plot_large %>% 
  pivot_longer(contains("cluster"), names_to = "cluster", values_to = "ancestry_prop")
qmat_plot_long$orig.pop <- factor(qmat_plot_long$orig.pop, levels = levels(data_sites$site))
qmat_plot_long <- 
  qmat_plot_long %>% 
  dplyr::arrange(orig.pop, Label, cluster, ancestry_prop) %>% 
  mutate(orig.pop = gsub('_', ' ', orig.pop)) %>% 
  mutate(orig.pop = gsub('W Australia', 'Western Australia', orig.pop)) %>% 
  mutate(orig.pop = factor(orig.pop, levels = LABELS)) 
  


### perso plot ----
gg_structure1 <- 
  ggplot(qmat_plot_long, aes(Label, ancestry_prop, fill = cluster)) +
  geom_col(linewidth = 0.01) + 
  facet_grid(~orig.pop, switch = "x", scales = "free", space = "free") +
  labs(x = "", y = "Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.15)) +
  scale_fill_viridis_d('Cluster', labels = c('Cluster 1', 'Cluster 2')) +
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    # axis.text.x = element_text(size=2, angle = 90, vjust = 0.5, hjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.background = element_rect(fill = 'white', color = 'white'),
    strip.background = element_rect(fill = 'white', color = "white")) 
gg_structure1





# ---- no Seychelles ----
## load ----
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "noSeychelles"

kmin <- 1
kmax <- 9
nrep <- 5
k <- 2
nloci <- 21905

dapc <- readRDS(paste0("intermediate/1_genetic_diversity/DAPC_popprior_output_", filters, "_", sites, "_", level, "_54PCA_8DA.RDS"))
structure <- readRDS(paste0("intermediate/1_genetic_diversity/STRUCTURE_", nloci, "loci_", sites, "_K", kmin, "-", kmax, "_r", nrep, ".RDS"))

scatter(dapc, scree.pca=TRUE)

##  scatter plot ----
# Plot ggplot - adapted from https://github.com/laurabenestan/TD_science_conservation

# setup axis
percent_explained <- dapc$eig / sum(dapc$eig) * 100
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("LD1 ({pretty_pe[1]}%)"),
            glue("LD2 ({pretty_pe[2]}%)"))


# get dapc values and metadata
dapc_geo <- 
  as.data.frame(dapc$ind.coord) %>% 
  rownames_to_column("id") %>% 
  dplyr::left_join(data_samples, by = "id") %>% 
  cbind(grp = dapc$grp)

# reorder labels by longitude
dapc_geo$grp <- 
  factor(dapc_geo$grp, 
         levels = unique(dapc_geo[order(dapc_geo$order),][[level]]))


# Plot
gg_dapc2 <- ggplot(dapc_geo, aes(x=LD1, y=LD2, color=grp)) +
  geom_vline(xintercept = 0, color = "darkgrey") +
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_point(size = 3, alpha = 0.5) +
  stat_ellipse(level = 0.67, show.legend = FALSE) +
  scale_color_manual(values = color_perso, labels = LABELS) +
  labs(x=labels[1], y=labels[2]) +
  theme_light() +
  theme(legend.position = "none") +
  labs(color = "")
gg_dapc2



## compolot ----

# create an object with membership probabilities
qmat <- gl.plot.structure(structure, K=k, colors_clusters = viridis(k))

### preparing data ----
qmat_plot_large <- 
  as.data.frame(qmat[[1]]) %>% 
  dplyr::select(c("orig.pop", "Label", contains("cluster"))) 
qmat_plot_long <- 
  qmat_plot_large %>% 
  pivot_longer(contains("cluster"), names_to = "cluster", values_to = "ancestry_prop")
qmat_plot_long$orig.pop <- factor(qmat_plot_long$orig.pop, levels = levels(data_sites$site))
qmat_plot_long <- 
  qmat_plot_long %>% 
  dplyr::arrange(orig.pop, Label, cluster, ancestry_prop) %>% 
  mutate(orig.pop = gsub('_', ' ', orig.pop)) %>% 
  mutate(orig.pop = gsub('W Australia', 'Western Australia', orig.pop)) %>% 
  mutate(orig.pop = factor(orig.pop, levels = LABELS)) 



### perso plot ----

# set up custom facet strips
facetstrips <-
  ggh4x::strip_nested(
    text_x = elem_list_text(size = c(8, 4), angle = 90, vjust = 0.5, hjust=1),
    by_layer_x = TRUE,
    clip = "off"
  )


# plot
gg_structure2 <- 
  ggplot(qmat_plot_long, aes(Label, ancestry_prop, fill = cluster)) +
  geom_col(linewidth = 0.01) + 
  facet_grid(~orig.pop, switch = "x", scales = "free", space = "free") +
  labs(x = "", y = "Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.15)) +
  scale_fill_viridis_d('Cluster', labels = c('Cluster 1', 'Cluster 2')) +
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    # axis.text.x = element_text(size=2, angle = 90, vjust = 0.5, hjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.background = element_rect(fill = 'white', color = 'white'),
    strip.background = element_rect(fill = 'white', color = "white")) 
gg_structure2




## {FIGURE 3} ####

gg <- 
  gg_dapc1 + gg_dapc2 + 
  gg_structure1 + gg_structure2 + 
  plot_layout(heights = c(4, 1), guides = 'collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')
gg

ggsave(paste0("results/1_genetic_diversity/_3_DAPC_STRUCTURE.png"),
       gg, 
       height = 8, width = 14, dpi = 500)
ggsave(paste0("results/1_genetic_diversity/_3_DAPC_STRUCTURE.pdf"),
       gg, 
       height = 8, width = 14)



