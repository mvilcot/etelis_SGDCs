
# ---- load ----
level = "site"

gd_alpha <- read.csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"))
gd_beta <- readRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))
BS <- readRDS(paste0("intermediate/1_genetic_diversity/basic_stats_", level, ".RDS"))


# ---- alpha GD by location ----

## Hs by loci ----
gd_alpha_loci <- 
  dplyr::as_tibble(BS[["Hs"]]) %>% 
  tidyr::pivot_longer(cols = everything(), 
               names_to = level, 
               values_to = "Hs")

# relevel sites
gd_alpha_loci[[level]] <- factor(gd_alpha_loci[[level]], levels = levels(data_samples[[level]]))

# plot
ggplot(gd_alpha_loci, aes(x=.data[[level]], y=.data[["Hs"]], fill = .data[[level]])) + 
  geom_boxplot() +
  scale_fill_manual(values = color_perso) +
  theme_light() +
  theme(legend.position = "none",
    axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

# save
ggsave(paste0("results/1_genetic_diversity/boxplot_alpha_Hs_", level, "_loci_hierfstat.png"),
       width = 8, height = 6)



## Hs by ind dartR ----
alpha_ind <- read_csv(paste0("results/1_genetic_diversity/gd_table_ind_dartR.csv"))

# add metadata
alpha_ind <-
  alpha_ind %>% 
  dplyr::rename(id = ind.name) %>% 
  full_join(data_samples, by= "id")

# plot
ggA <- 
  ggplot(alpha_ind, aes(x=.data[[level]], y=.data[["Ho"]], fill = .data[[level]])) + 
  geom_boxplot() +
  scale_fill_manual(values = color_perso) +
  theme_light() +
  theme(legend.position = "none",
    axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave(paste0("results/1_genetic_diversity/boxplot_alpha_Ho_", level, "_ind_dartR.png"),
       ggA,
       width = 8, height = 6)


# ---- Hs ~ Fst pop specific ----
# Pearson, lm
corP <- cor.test(gd_alpha$Hs, gd_alpha$popFst.WG)
# model <- summary(lm(Hs ~ popFst.WG, data = gd_alpha))

# relevel sites
gd_alpha[[level]] <- factor(gd_alpha[[level]], levels = levels(data_samples[[level]]))

# plot
gd_alpha <- 
  gd_alpha %>% 
  mutate(Labels = LABELS) 
ggplot(gd_alpha, aes(x=popFst.WG, y=Hs, color = .data[[level]])) + 
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(data = gd_alpha[gd_alpha$site %in% c("Seychelles","Christmas_Island","Hawaii","W_Australia"),],
                  aes(label=.data[["Labels"]]), box.padding = 0.5, size = 3,
                  bg.color = "grey70", bg.r = 0.02, vjust = 0.5, hjust = 0.5, point.padding = 0.5, force = 100,
                  show.legend  = FALSE) +
  scale_color_manual('', values = color_perso, 
                     labels = LABELS) +
  annotate('text', x=max(gd_alpha$popFst.WG), y=max(gd_alpha$Hs), 
           hjust=1, vjust=0.5, size=3,
           label=paste0("r Pearson = ", signif(corP$estimate, 3), "\n p = ", signif(corP$p.value, 1))) +
  xlab('site-specific Fst') +
  theme_light()

## {FIGURE 2} ####
ggsave(paste0("results/1_genetic_diversity/_2_plot_Hs_Fst_pop_specific_", level, "_light.png"),
       width = 6, height = 4, dpi = 500)
ggsave(paste0("results/1_genetic_diversity/_2_plot_Hs_Fst_pop_specific_", level, "_light.pdf"),
       width = 6, height = 4)



# Hs ~ longitude ----
## load
gd_alpha$site <- factor(gd_alpha$site, ordered = T)

df <-
  data_sites %>% 
  left_join(gd_alpha)

## shift coordinates
df <- shift.lon(df)# Pearson, lm

## correlation
corP <- cor.test(df$longitude, df$Ho)

ggplot(df, aes(x=longitude, y=Ho, color = .data[[level]])) + 
  # geom_smooth(method='lm', formula=y~x, colour = "grey35") +
  geom_point() +
  scale_color_manual(values = color_perso) +
  theme_light()+
  # geom_text_repel(aes(label=.data[[level]]), ) +
  annotate('text', x=max(df$longitude), y=max(df$Ho), 
           hjust=1, vjust=1, size=4.5,
           label=paste0("Pearson = ", round(corP$estimate, 4),
                        "\n p = ", round(corP$p.value, 4)))


# save
ggsave(paste0("results/1_genetic_diversity/_alpha_gd_Ho_longitude_site.png"),
       width = 8, height = 6)




# ---- comparison beta GD metrics ----
list_GDbeta <- list()

for (metricGD in names(gd_beta)){
  
  # get distance matrix
  mat_GDbeta <- as.matrix(gd_beta[[metricGD]])
  
  # order rows alphabetically
  mat_GDbeta <- mat_GDbeta[order(rownames(mat_GDbeta)), order(colnames(mat_GDbeta))]
  
  # pivot longer distance matrix
  melt_GDbeta <- melt.dist(dist = mat_GDbeta, metric = metricGD)
  
  # put into list
  list_GDbeta[[metricGD]] <- melt_GDbeta
  
}

# merge two distance matrix into one df
merge_gd_beta <- 
  plyr::join_all(list_GDbeta,
                 by = c(paste0(level, "1"), paste0(level, "2")))

merge_gd_beta$site <- paste(merge_gd_beta$site1,
                            merge_gd_beta$site2,
                            sep = "-")

merge_gd_beta$Christmas <- grepl("Christmas", merge_gd_beta$site)


ggplot(merge_gd_beta, aes(Fst, GstPP.hed, color = Christmas)) +
  geom_point() 

ggplot(merge_gd_beta, aes(D.Jost, GstPP.hed, color = Christmas)) +
  geom_point()

ggplot(merge_gd_beta, aes(Fst, D.Jost, color = Christmas)) +
  geom_point()


# library(plotly)
# plot_ly(x=merge_gd_beta[["GstPP.hed"]], 
#         y=merge_gd_beta[["JostD"]], 
#         z=merge_gd_beta[["Fst"]], 
#         type="scatter3d", mode="markers", size = 0.5)
# 
#  
# 
# library(corrplot)
# temp <- cor(merge_gd_beta[, names(gd_beta)])
# corrplot(temp)
# 


# Hs ~ dist to CT ----
## !!!!!!!! TO DO !!!!!!!!!!!!!!!! ----



