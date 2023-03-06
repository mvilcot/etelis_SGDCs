
## ---- load ----
level = "site"

gd_beta <- readRDS(paste0("results/01_genetic_diversity/gd_list_pairwise_", level, ".RDS"))
bs <- readRDS("intermediate/01_genetic_diversity/basic_stats.RDS")



## ---- alpha GD by site ----

gd_alpha_loci <- 
  as_tibble(bs[["Hs"]]) %>% 
  pivot_longer(cols = everything(), 
               names_to = level, 
               values_to = "Hs")

ggplot(gd_alpha_loci, aes(x=site, y=.data[["Hs"]])) + 
  geom_boxplot()



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

ggplot(merge_gd_beta, aes(.data[["Fst"]], .data[["GstPP.hed"]])) +
  geom_point()

ggplot(merge_gd_beta, aes(.data[["JostD"]], .data[["GstPP.hed"]])) +
  geom_point()

ggplot(merge_gd_beta, aes(.data[["Fst"]], .data[["GstPP.hed"]])) +
  geom_point()


library(plotly)
plot_ly(x=merge_gd_beta[["GstPP.hed"]], 
        y=merge_gd_beta[["JostD"]], 
        z=merge_gd_beta[["Fst"]], 
        type="scatter3d", mode="markers", size = 0.5)



library(corrplot)
temp <- cor(merge_gd_beta[, names(gd_beta)])
corrplot(temp)

