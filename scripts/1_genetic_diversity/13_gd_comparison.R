
# ---- load ----
level = "site"

gd_alpha <- read.csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"))
gd_beta <- readRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))
BS <- readRDS(paste0("intermediate/1_genetic_diversity/basic_stats_", level, ".RDS"))


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



