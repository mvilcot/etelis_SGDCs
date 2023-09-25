
# function to Melt distance matrix ----
melt.dist <- function(distmat, colname) {
  if(class(distmat)[1] == "dist") {distmat <- as.matrix(distmat)}
  distmat[upper.tri(distmat, diag = T)] <- NA
  distmat <- 
    as.data.frame(distmat) %>% 
    rownames_to_column("site1") %>% 
    pivot_longer(cols = -"site1", 
                 names_to = "site2", 
                 values_to = colname) %>% 
    na.omit()
  
  return(distmat)
}


# ---- melt ----
dist_melt <- list()
for (i in 1:length(dist_mat)){
  dist_melt[[i]] <- melt.dist(distmat = dist_mat[[i]], metric = names(dist_mat)[[i]])
}

# ---- merge into a table ----
# full join
dist_merge <- 
  dist_melt %>%
  reduce(full_join, by = c("site1", "site2")) 

# create one column for both locations info
dist_merge$site <- 
  paste(dist_merge[["site1"]], 
        dist_merge[["site2"]],
        sep = "-")

# relocate to first column
dist_merge <- 
  dist_merge %>% 
  relocate("site")





# load data ----
dist_merge <-
  read_csv(paste0("results/3_distance_metrics/dist_geo_envt_res17-4_gd_sd_", comm_delin, ".csv"))



# linear model ----
VARIABLE <- dist_merge[["Fst"]]
FACTOR <- dist_merge[,c("environment", "seadist")]

VARIABLE <- decostand(VARIABLE, method="standardize")[,1]
FACTOR <- decostand(FACTOR, method="standardize")

LM.model <- lm(VARIABLE~., data=FACTOR)
summary(LM.model)


# MRM ----
MRM(VARIABLE~., data=FACTOR)

# plot 
dwplot(LM.model)




# anova table ----
ANOVA.model <- anova(LM.model)
df <- data.frame(variable = rownames(ANOVA.model),
                 variance = ANOVA.model$`Sum Sq` / sum(ANOVA.model$`Sum Sq`),
                 x=NA)

# base plot ----
barplot(df$variance, horiz = TRUE, col=c("black", "red", "darkgrey"),
        names.arg = df$variable,
        main = "variance decomposition", xlab = "% of the model")

# ggplot ----
ggplot(df, aes(fill=variable, y=variance, x=x)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_flip()


