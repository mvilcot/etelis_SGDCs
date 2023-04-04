

## ---- load ----
level = "site"

gd_alpha <- read.csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"))
gd_beta <- readRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))
BS <- readRDS(paste0("intermediate/1_genetic_diversity/basic_stats_", level, ".RDS"))

# order_sites <- 
#   data_sites[[level]] %>% 
#   levels(unique(data_sites[order(data_sites$order),][[level]])) %>% 
#   droplevels()
# data_samples$site




## ---- alpha GD by location ----

# get Hs by loci by site
gd_alpha_loci <- 
  as_tibble(BS[["Hs"]]) %>% 
  pivot_longer(cols = everything(), 
               names_to = level, 
               values_to = "Hs")

# relevel sites
gd_alpha_loci[[level]] <- factor(gd_alpha_loci[[level]], levels = levels(data_samples[[level]]))

# plot
ggplot(gd_alpha_loci, aes(x=.data[[level]], y=.data[["Hs"]], fill = .data[[level]])) + 
  geom_boxplot() +
  scale_fill_manual(values = color_perso) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# save
ggsave(paste0("results/1_genetic_diversity/plot_alpha_Hs_by_", level, ".png"),
       width = 10, height = 8)




## ---- Hs ~ Fst pop specific ----

library(ggrepel)

# Pearson, lm
corP <- cor(gd_alpha$Hs, gd_alpha$popFst.WG)
model <- summary(lm(Hs ~ popFst.WG, data = gd_alpha))

# relevel sites
gd_alpha[[level]] <- factor(gd_alpha[[level]], levels = levels(data_samples[[level]]))

# plot
ggplot(gd_alpha, aes(x=popFst.WG, y=Hs, color = .data[[level]])) + 
  geom_smooth(method='lm', formula=y~x, colour = "grey35") +
  geom_point() +
  scale_color_manual(values = color_perso) +
  geom_text_repel(aes(label=.data[[level]]), ) +
  annotate('text', x=max(gd_alpha$popFst.WG), y=max(gd_alpha$Hs), 
           hjust=1, vjust=1, size=4.5,
           label=paste0("LM adjusted R² = ", round(model$adj.r.squared, 3), 
                        ", P = ", round(coef(model)[2,4], 4),
                        "\n Pearson = ", round(corP, 4)))

# save
ggsave(paste0("results/1_genetic_diversity/plot_Hs_Fst_pop_specific_", level, ".png"),
       width = 8, height = 6)





## ---- beta GD metrics ----
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
