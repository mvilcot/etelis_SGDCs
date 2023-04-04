
## ---- load data ----

# parameters
level = "site"
comm_delin = "taxonomic_scale"

# read GD and SD genetic diversity
gd_beta <- readRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))
sd_beta <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_", level, "_", comm_delin, ".RDS"))

# communities delineation
list_communities <- readRDS(paste0("intermediate/2_species_diversity/List_community_", comm_delin, ".RDS"))

# parameters
metricSD = "beta.jtu"
metricGD = "Fst"
comm = names(list_communities)[1]



## ---- Matrix distance ----
mat_lcdist <- 
  read.csv(file = "intermediate/3_distance_decay/least_cost_distance.csv", row.names = 1)


## ---- Matrix distance ----
mat_IBR <- 
  read.csv("intermediate/5_re_Lesturgie/least_cost_distance_IBR_23_100.csv", row.names = 1)



## ---- Matrix environmental variables ----
# read variables 
site_variables <- read.csv("intermediate/3_distance_decay/bio_oracle_variables.csv",
                           row.names = 1)

env <- site_variables[, c(1, 10)] 
pca_env <- dudi.pca(df = env, scannf = FALSE, nf = 3) 

mat_env <- 
  vegdist(pca_env$li, method = "euclidean", na.rm = TRUE) %>% # euclidean distance based on the 3 PCA axes
  as.matrix()

pca_coord <- 
  pca_env$li %>% 
  rownames_to_column(level)

ggplot(pca_coord, aes(Axis1, Axis2, label = .data[[level]])) +
  geom_point() +
  geom_text(hjust=0, vjust=0)


# ## ---- merge to GD and SD ----
# merge_beta <- read.csv("intermediate/3_distance_decay/temp_merge_beta.csv")
# 
# # pivot longer distance matrix
# melt_env <- 
#   melt.dist(dist = mat_env, metric = "distenv") 
# 
# # merge two distance matrix into one df
# merge_beta <- 
#   merge_beta %>% 
#   left_join(melt_env,
#             by = c(paste0(level, "1"), paste0(level, "2"))) 
# 
# 
# ## ---- SGDCs decomposition ----
# SGDC.decomp(SD = merge_beta[[metricSD]], 
#             GD = merge_beta[[metricGD]], 
#             FACTOR = merge_beta[,c("lcdist","distenv")])




## ---- setup beta data ----
# get distance matrix
mat_SDbeta <- as.matrix(sd_beta[[comm]][[metricSD]])
mat_GDbeta <- as.matrix(gd_beta[[metricGD]])

# keep only sampling sites present in both GD and SD
mat_SDbeta <- mat_SDbeta[rownames(mat_SDbeta) %in% rownames(mat_GDbeta),
                         colnames(mat_SDbeta) %in% colnames(mat_GDbeta)]
mat_GDbeta <- mat_GDbeta[rownames(mat_GDbeta) %in% rownames(mat_SDbeta),
                         colnames(mat_GDbeta) %in% colnames(mat_SDbeta)]
mat_env <- mat_env[rownames(mat_env) %in% rownames(mat_SDbeta),
                   colnames(mat_env) %in% colnames(mat_SDbeta)]
mat_lcdist <- mat_lcdist[rownames(mat_lcdist) %in% rownames(mat_SDbeta),
                         colnames(mat_lcdist) %in% colnames(mat_SDbeta)]
mat_IBR <- mat_IBR[rownames(mat_IBR) %in% rownames(mat_SDbeta),
                   colnames(mat_IBR) %in% colnames(mat_SDbeta)]

# order rows alphabetically
mat_SDbeta <- as.dist(mat_SDbeta[order(rownames(mat_SDbeta)), order(colnames(mat_SDbeta))])
mat_GDbeta <- as.dist(mat_GDbeta[order(rownames(mat_GDbeta)), order(colnames(mat_GDbeta))])
mat_env <- as.dist(mat_env[order(rownames(mat_env)), order(colnames(mat_env))])
mat_lcdist <- as.dist(mat_lcdist[order(rownames(mat_lcdist)), order(colnames(mat_lcdist))])
mat_IBR <- as.dist(mat_IBR[order(rownames(mat_IBR)), order(colnames(mat_IBR))])






## ---- MRM ----

MRM(mat_SDbeta ~ mat_lcdist + mat_env , nperm = 9999)
MRM(mat_SDbeta ~ mat_IBR + mat_lcdist + mat_env , nperm = 9999)

MRM(mat_GDbeta ~ mat_lcdist + mat_env , nperm = 9999)
MRM(mat_GDbeta ~ mat_IBR + mat_env , nperm = 9999)
MRM(mat_GDbeta ~ mat_IBR)
MRM(mat_GDbeta ~ mat_lcdist)


cor(mat_lcdist, mat_env)
cor(mat_lcdist, mat_IBR)
plot(mat_lcdist, mat_IBR)


# 
# ## ----------------------- Robuchon GD SD decomposition ---------------------- 
# 
# 
# ## ---- step 1: perform the lm and compute standardised partial regression coefficients ----
# lm.species <- lm(list.Betasim_fish_species[[i]] ~ list.basinsdist_species[[i]]
#                  + list.env_species[[i]] +  list.area_species[[i]])
# lm.table_betaSD$Intercept[i] <-  lm.species$coefficients[1]
# lm.table_betaSD$Geo[i] <-  lm.species$coefficients[2]*sd(list.basinsdist_species[[i]])/sd(list.Betasim_fish_species[[i]])
# lm.table_betaSD$Env[i] <-  lm.species$coefficients[3]*sd(list.env_species[[i]])/sd(list.Betasim_fish_species[[i]])
# lm.table_betaSD$Area[i] <-  lm.species$coefficients[4]*sd(list.area_species[[i]])/sd(list.Betasim_fish_species[[i]])
# lm.table_betaSD$R2[i] <-  summary(lm.species)$r.squared
# lm.table_betaSD$pval.R2[i] <- pf(summary(lm.species)$fstatistic[1], summary(lm.species)$fstatistic[2], 
#                                  summary(lm.species)$fstatistic[3],lower.tail=FALSE) # OK until here
# res_betaSD[[i]] <- lm.species$residuals # save residuals
# 
# ## ---- step 2: permute row and columns of the distance matrix 1999 times ---- 
# # run linear regression for each of these 1999 permuted distance matrices
# # and compute null distributions of pseudo-t associated to each regression coefficients
# perm.Betasim <- replicate(1999, permute(list.Betasim_fish_species[[i]], sample(1:attributes(list.Betasim_fish_species[[i]])$Size)))
# null.distr <- data.frame(row.names = 1:1999)
# 
# for (j in 1:1999)
# {
#   perm.lm <- lm(perm.Betasim[,j] ~ as.numeric(list.basinsdist_species[[i]]) + as.numeric(list.env_species[[i]]) + 
#                   as.numeric(list.area_species[[i]]))
#   # calculate pseudo-t
#   tcorr <- sqrt(1 - summary(perm.lm)$r.squared)
#   # compute null distrib
#   null.distr$Intercept[j] <- perm.lm$coefficients[1] / tcorr
#   null.distr$Geo[j] <- perm.lm$coefficients[2] / tcorr
#   null.distr$Env[j] <- perm.lm$coefficients[3] / tcorr
#   null.distr$Area[j] <- perm.lm$coefficients[4] / tcorr
# }
# # include unpermuted data in the null distribution
# unpermuted <- c(lm.species$coefficients[1]/tcorr, lm.species$coefficients[2]/tcorr, lm.species$coefficients[3]/tcorr, 
#                 lm.species$coefficients[4]/tcorr)
# null.distr <- rbind(unpermuted, null.distr)
# 
# # ---- step 3: compare obs coef to null distributions to compute p-values ----
# lm.table_betaSD$pval.Intercept[i] <- length(null.distr$Intercept[abs(null.distr$Intercept) >= abs(null.distr$Intercept[1])])/2000 # 2-sided
# lm.table_betaSD$pval.Geo[i] <- length(null.distr$Geo[null.distr$Geo > null.distr$Geo[1]])/2000 # 1-sided, positive effect
# lm.table_betaSD$pval.Env[i] <- length(null.distr$Env[null.distr$Env > null.distr$Env[1]])/2000 # 1-sided, positive effect
# lm.table_betaSD$pval.Area[i] <- length(null.distr$Area[null.distr$Area < null.distr$Area[1]])/2000 # 1-sided, negative effect
# 
# 
# 
# 
# 
# 
# ## ---- Partial effects on SD ----
# # distributions
# SD.effects <- data.frame (lm.table_betaSD$Geo, lm.table_betaSD$Env, lm.table_betaSD$Area)
# colnames(SD.effects) <- c("Geographic distance", "Environmental distance", "Harmonic mean area")
# SD.effects <- cbind(species = rownames(lm.table_betaSD), SD.effects)
# mcomb <- melt (SD.effects, id.vars = "species", 
#                measure.vars = c("Geographic distance", "Environmental distance", "Harmonic mean area"))
# g.SD <- ggplot(mcomb, aes(variable, value, fill = variable))
# g.SD <- g.SD + geom_hline(yintercept = 0, linetype ="dashed") + geom_violin(scale = "area") +
#   theme_bw(base_size = 8) +
#   # geom_boxplot(width = 0.1) +
#   stat_summary(fun.data = data_summary, geom = "crossbar", width = 0.1, color = "black", fill = "white") +
#   labs(title = "(a) Multiple regressions on taxonomic differentiation", x = "", y ="standardised regression coefficient") +
#   scale_fill_brewer(palette = "Dark2", guide = FALSE) +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize = 1, alpha = 0.5, binwidth = 0.1, fill = "black") +
#   theme(plot.title = element_text(size = 10), axis.text.x = element_text(color = brewer.pal(3, name = "Dark2")))
# g.SD
# 
# # tests
# SDgeo_pvalmean <- mean(lm.table_betaSD$pval.Geo)
# SDgeo_null.distrib <- c()
# for (i in 1:10000)
# {
#   SDgeo_null.distrib[i] <- mean(runif(nrow(lm.table_betaSD)))
# }
# ptest_SDgeo <- length(SDgeo_null.distrib[which(SDgeo_null.distrib <= SDgeo_pvalmean)])/10000
# 
# SDenv_pvalmean <- mean(lm.table_betaSD$pval.Env)
# SDenv_null.distrib <- c()
# for (i in 1:10000)
# {
#   SDenv_null.distrib[i] <- mean(runif(nrow(lm.table_betaSD)))
# }
# ptest_SDenv <- length(SDenv_null.distrib[which(SDenv_null.distrib <= SDenv_pvalmean)])/10000
# 
# SDarea_pvalmean <- mean(lm.table_betaSD$pval.Area)
# SDarea_null.distrib <- c()
# for (i in 1:10000)
# {
#   SDarea_null.distrib[i] <- mean(runif(nrow(lm.table_betaSD)))
# }
# ptest_SDarea <- length(SDarea_null.distrib[which(SDarea_null.distrib <= SDarea_pvalmean)])/10000
# 
# df_null.distrib <- data.frame(SDgeo_null.distrib, SDenv_null.distrib, SDarea_null.distrib)
# colnames(df_null.distrib) <- c("Geo", "Env", "Area")
# 
# t.SDgeo <- ggplot(df_null.distrib, aes(Geo)) + geom_histogram(binwidth = 0.02, color = "black", fill = "grey", size = 0.1) +
#   theme_bw(base_size = 8) + geom_vline(xintercept = SDgeo_pvalmean, linetype ="dashed", size = 0.5) +
#   scale_x_continuous(limits = c(0.1, 0.9)) + scale_y_continuous(limits = c(0, 1500)) +
#   labs(x = "null distribution of p-values", y ="", title = "(a) Multiple regressions on taxonomic differentiation",
#        subtitle = "Geographic distance") +
#   theme(plot.subtitle = element_text(size = 8, color = "#1B9E77"),
#         plot.margin = margin(0.05, 0.05, 0, 0, "cm"), axis.text.y = element_blank(),
#         axis.ticks = element_blank()) +
#   annotate(geom = "text", x = 0.6, y = 1500, size = 2,
#            label = paste0("mean observed p-value = ", round(SDgeo_pvalmean, 4), ", p-value(test) = ", round(ptest_SDgeo, 4))) 
# t.SDgeo
# 
# t.SDenv <- ggplot(df_null.distrib, aes(Env)) + geom_histogram(binwidth = 0.02, color = "black", fill = "grey", size = 0.1) +
#   theme_bw(base_size = 8) + geom_vline(xintercept = SDenv_pvalmean, linetype ="dashed", size = 0.5) +
#   scale_x_continuous(limits = c(0.1, 0.9)) + scale_y_continuous(limits = c(0, 1500)) +
#   labs(x = "null distribution of p-values", y ="", subtitle = "Environmental distance") +
#   theme(plot.subtitle = element_text(size = 8, color = "#D95F02"),
#         plot.margin = margin(0.05, 0.05, 0, 0, "cm"), axis.text.y = element_blank(),
#         axis.ticks = element_blank()) +
#   annotate(geom = "text", x = 0.6, y = 1500, size = 2,
#            label = paste0("mean observed p-value = ", round(SDenv_pvalmean, 4), ", p-value(test) = ", round(ptest_SDenv, 4))) 
# t.SDenv
# 
# t.SDarea <- ggplot(df_null.distrib, aes(Area)) + geom_histogram(binwidth = 0.02, color = "black", fill = "grey", size = 0.1) +
#   theme_bw(base_size = 8) + geom_vline(xintercept = SDarea_pvalmean, linetype ="dashed", size = 0.5) +
#   scale_x_continuous(limits = c(0.1, 0.9)) + scale_y_continuous(limits = c(0, 1500)) +
#   labs(x = "null distribution of p-values", y ="", 
#        subtitle = "Harmonic mean area") +
#   theme(plot.subtitle = element_text(size = 8, color = "#7570B3"),
#         plot.margin = margin(0.05, 0.05, 0, 0, "cm"), axis.text.y = element_blank(),
#         axis.ticks = element_blank()) +
#   annotate(geom = "text", x = 0.6, y = 1500, size = 2,
#            label = paste0("mean observed p-value = ", round(SDarea_pvalmean, 4), ", p-value(test) = ", round(ptest_SDarea, 4))) 
# t.SDarea
# 
# ## ---- Partial effects of Geo, Env & Area on GD ----
# # distributions
# GD.effects <- data.frame (lm.table_betaGD$Geo, lm.table_betaGD$Env, lm.table_betaGD$Area)
# colnames(GD.effects) <- c("Geographic distance", "Environmental distance", "Harmonic mean area")
# GD.effects <- cbind(species = rownames(lm.table_betaGD), GD.effects)
# mcomb <- melt (GD.effects, id.vars = "species", 
#                measure.vars = c("Geographic distance", "Environmental distance", "Harmonic mean area"))
# g.GD <- ggplot(mcomb, aes(variable, value, fill = variable))
# g.GD <- g.GD + geom_hline(yintercept = 0, linetype ="dashed") + geom_violin(scale = "area") +
#   theme_bw(base_size = 8) +
#   # geom_boxplot(width = 0.1) +
#   stat_summary(fun.data = data_summary, geom = "crossbar", width = 0.1, color = "black", fill = "white") +
#   labs(title = "(b) Multiple regressions on genetic differentiation", x = "", y ="standardised regression coefficient") +
#   scale_fill_brewer(palette = "Dark2", guide = FALSE) +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize = 1, alpha = 0.5, binwidth = 0.1, fill = "black") +
#   theme(plot.title = element_text(size = 10), axis.text.x = element_text(color = brewer.pal(3, name = "Dark2")))
# g.GD
# 
# # tests with Hedrick
# GDgeo_pvalmean <- mean(lm.table_betaGD$pval.Geo)
# GDgeo_null.distrib <- c()
# for (i in 1:10000)
# {
#   GDgeo_null.distrib[i] <- mean(runif(nrow(lm.table_betaGD)))
# }
# ptest_GDgeo <- length(GDgeo_null.distrib[which(GDgeo_null.distrib <= GDgeo_pvalmean)])/10000
# 
# GDenv_pvalmean <- mean(lm.table_betaGD$pval.Env)
# GDenv_null.distrib <- c()
# for (i in 1:10000)
# {
#   GDenv_null.distrib[i] <- mean(runif(nrow(lm.table_betaGD)))
# }
# ptest_GDenv <- length(GDenv_null.distrib[which(GDenv_null.distrib <= GDenv_pvalmean)])/10000
# 
# GDarea_pvalmean <- mean(lm.table_betaGD$pval.Area)
# GDarea_null.distrib <- c()
# for (i in 1:10000)
# {
#   GDarea_null.distrib[i] <- mean(runif(nrow(lm.table_betaGD)))
# }
# ptest_GDarea <- length(GDarea_null.distrib[which(GDarea_null.distrib <= GDarea_pvalmean)])/10000
# 
# df_null.distrib <- data.frame(GDgeo_null.distrib, GDenv_null.distrib, GDarea_null.distrib)
# colnames(df_null.distrib) <- c("Geo", "Env", "Area")
# 
# t.GDgeo <- ggplot(df_null.distrib, aes(Geo)) + geom_histogram(binwidth = 0.02, color = "black", fill = "grey", size = 0.1) +
#   theme_bw(base_size = 8) + geom_vline(xintercept = GDgeo_pvalmean, linetype ="dashed", size = 0.5) +
#   scale_x_continuous(limits = c(0.1, 0.9)) + scale_y_continuous(limits = c(0, 1500)) +
#   labs(x = "null distribution of p-values", y ="", title = "(b) Multiple regressions on genetic differentiation (G''ST)",
#        subtitle = "Geographic distance") +
#   theme(plot.subtitle = element_text(size = 8, color = "#1B9E77"),
#         plot.margin = margin(0.05, 0.05, 0, 0, "cm"), axis.text.y = element_blank(),
#         axis.ticks = element_blank()) +
#   annotate(geom = "text", x = 0.6, y = 1500, size = 2,
#            label = paste0("mean observed p-value = ", round(GDgeo_pvalmean, 4), ", p-value(test) = ", round(ptest_GDgeo, 4))) 
# t.GDgeo
# 
# t.GDenv <- ggplot(df_null.distrib, aes(Env)) + geom_histogram(binwidth = 0.02, color = "black", fill = "grey", size = 0.1) +
#   theme_bw(base_size = 8) + geom_vline(xintercept = GDenv_pvalmean, linetype ="dashed", size = 0.5) +
#   scale_x_continuous(limits = c(0.1, 0.9)) + scale_y_continuous(limits = c(0, 1500)) +
#   labs(x = "null distribution of p-values", y ="", subtitle = "Environmental distance") +
#   theme(plot.subtitle = element_text(size = 8, color = "#D95F02"),
#         plot.margin = margin(0.05, 0.05, 0, 0, "cm"), axis.text.y = element_blank(),
#         axis.ticks = element_blank()) +
#   annotate(geom = "text", x = 0.6, y = 1500, size = 2,
#            label = paste0("mean observed p-value = ", round(GDenv_pvalmean, 4), ", p-value(test) = ", round(ptest_GDenv, 4))) 
# t.GDenv
# 
# t.GDarea <- ggplot(df_null.distrib, aes(Area)) + geom_histogram(binwidth = 0.02, color = "black", fill = "grey", size = 0.1) +
#   theme_bw(base_size = 8) + geom_vline(xintercept = GDarea_pvalmean, linetype ="dashed", size = 0.5) +
#   scale_x_continuous(limits = c(0.1, 0.9)) + scale_y_continuous(limits = c(0, 1500)) +
#   labs(x = "null distribution of p-values", y ="", 
#        subtitle = "Harmonic mean area") +
#   theme(plot.subtitle = element_text(size = 8, color = "#7570B3"),
#         plot.margin = margin(0.05, 0.05, 0, 0, "cm"), axis.text.y = element_blank(),
#         axis.ticks = element_blank()) +
#   annotate(geom = "text", x = 0.6, y = 1500, size = 2,
#            label = paste0("mean observed p-value = ", round(GDarea_pvalmean, 4), ", p-value(test) = ", round(ptest_GDarea, 4))) 
# t.GDarea
# 
# 
