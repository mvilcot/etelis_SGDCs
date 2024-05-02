
# ---- Read SNPs data set ----
# parameters
## :::: To adapt if allsites or noSeychelles ####
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "allsites" 

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL, #"Seychelles"
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)


# ---- Convert to LEA data format ----

# SNP presence/absence lfmm (package LEA, SilicoDArT)
dartR::gl2geno(genlight,
        outfile = paste0("GenoLEA_Etelis_coruscans_ordered_", sites),
        outpath='intermediate/1_genetic_diversity')



# ---- STRUCTURE ----
# install STUCTURE (non GUI version) from here: https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html

## setup genlight ----
genlightLL <- genlight

# add lon/lat (when location by site)
genlightLL@other$latlon <-
  genlightLL@other$ind.metrics %>%
  left_join(shift.lon(data_sites), by = "site") %>%
  dplyr::select(latitude, longitude) %>%
  dplyr::rename(lat = latitude) %>%
  dplyr::rename(lon = longitude)


## map lat lon ----
gl.map.interactive(genlightLL)

# change levels to alphabetic, otherwise issue with gl.map.structure
genlightLL$pop <- factor(as.character(genlightLL$pop))


## run structure ----
kmin <- 1
kmax <- 9
nrep <- 5
nloci <- genlightLL$n.loc

structure <- 
  dartR::gl.run.structure(
    genlightLL,
    k.range = kmin:kmax,
    num.k.rep = nrep,
    # pop.prior = "usepopinfo",
    exec = "C:/Users/vilcot/Documents/Structure2.3.4/console/structure.exe",
    plot.out = FALSE)

structure %>% saveRDS(paste0("intermediate/1_genetic_diversity/STRUCTURE_", nloci, "loci_", sites, "_K", kmin, "-", kmax, "_r", nrep, ".RDS"))


## output ---- 
structure <- readRDS(paste0("intermediate/1_genetic_diversity/STRUCTURE_", nloci, "loci_", sites, "_K", kmin, "-", kmax, "_r", nrep, ".RDS"))


# evanno plot
ev <- gl.evanno(structure)

evgg <- # arrange plots
  wrap_plots(ev$plots) +
  patchwork::plot_annotation(tag_levels = 'a', # add tags
                             tag_prefix = '(',
                             tag_suffix = ')',
                             tag_sep = "")

evgg <- # add common x axis
  patchwork::wrap_elements(panel = evgg) +
  labs(tag = "K") + 
  theme(plot.tag = element_text(size = rel(1)),
        plot.tag.position = "bottom")

evgg
ggsave(paste0("results/1_genetic_diversity/STRUCTURE_evanno_", nloci, "loci_", sites, "_K", kmin, "-", kmax, "_r", nrep, ".png"),
       evgg, 
       height = 6, width = 8)


# automatic ancestry plot
k <- 2

# proba matrix
qmat <- gl.plot.structure(structure, K=k, colors_clusters = viridis(k))

# map barplot
gl.map.structure(qmat, genlightLL, K=k, scalex=0.5, scaley=0.2) 



## plot structure ----

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
  dplyr::arrange(orig.pop, ancestry_prop, cluster)



###  plot ----
gg <- 
  ggplot(qmat_plot_long, aes(Label, ancestry_prop, fill = cluster)) +
  geom_col(width = 1.1) +
  facet_grid(~orig.pop, switch = "x", scales = "free", space = "free") +
  labs(x = "", y = "Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.15)) +
  scale_fill_viridis_d('', labels = c('Cluster 1', 'Cluster 2')) +
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.background = element_rect(fill = 'white', color = 'white'),
    strip.background = element_rect(fill = 'white', color = "white"))

gg




### map ----
library(scatterpie)

# Proportion of ancestry in each site
prop_cluster <-
  qmat_plot_large %>% 
  dplyr::rename(site = orig.pop) %>% 
  group_by(site) %>% 
  dplyr::summarise(cluster1 = sum(cluster1), 
                   cluster2 = sum(cluster2), 
                   # cluster3 = sum(cluster3), 
                   # cluster4 = sum(cluster4),
                   # cluster5 = sum(cluster5),
                   .groups = "keep") 
prop_cluster$site <- ordered(prop_cluster$site, levels=levels(data_sites$site))

prop_cluster <-
  prop_cluster %>% 
  left_join(data_sites, by = "site")

# load maps --
# full map
map1 <-
  fortify(maps::map(fill=TRUE, plot=FALSE)) %>% 
  as_tibble()

# add map +360
map2 <- 
  map1 %>% 
  mutate(long = long + 360,
         group = group + max(map1$group) + 1)

# crop lon & lat extent
map <- 
  rbind(map1, map2) %>% 
  filter(long > 20  & long < 210 & lat <48 & lat > -48)

gg2 <- 
  ggplot() +
  ## Countries
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
  ## Sites
  geom_scatterpie(data = shift.lon(prop_cluster), 
             aes(x = longitude, y = latitude), #r = number_samples/10
             cols = c("cluster1", "cluster2"), # , "cluster3", "cluster4", "cluster5"
             pie_scale = 1.5, color=NA) + 
  ggrepel::geom_text_repel(data = shift.lon(prop_cluster), show.legend = FALSE,
                           aes(x = longitude, y = latitude, color = site, label = site),
                           hjust=0.5, vjust=0, max.overlaps = 5, nudge_x = -10,
                           bg.color = "grey70", bg.r = 0.02) +
  ## Theme
  scale_color_viridis(discrete = T) +
  scale_fill_viridis('', discrete = T, direction = 1, labels = c('Cluster 1', 'Cluster 2')) +
  theme_minimal() +
  theme(plot.background = element_rect(fill="white", color = "white")) +
  labs(x = "", y = "") +
  coord_equal()
  
gg1


#### {FIGURE S3} ####
## save 
gg12 <- 
  gg1 / gg2 +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")") 
ggsave(paste0("results/1_genetic_diversity/_S3_STRUCTURE_mappiechart_K", k, "_propancestry2.png"),
       gg12,
       height = 10, width = 10, dpi = 500)



