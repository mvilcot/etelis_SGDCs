## ---- read variables ----
site_variables <- read.csv("intermediate/03_distance_decay/bio_oracle_variables.csv",
                           row.names = 1)

## ---- PCA on the environmental variables ----
env <- site_variables[, c(1, 10)] 
pca.env <- dudi.pca(df = env, scannf = FALSE, nf = 3) 

mat_env <- 
  vegdist(pca.env$li, method = "euclidean", na.rm = TRUE) %>% # euclidean distance based on the 3 PCA axes
  as.matrix()

## ---- merge to GD and SD ----
merge_beta <- read.csv("intermediate/03_distance_decay/temp_merge_beta.csv")

# pivot longer distance matrix
melt_env <- 
  melt.dist(dist = mat_env, metric = "distenv") 

# merge two distance matrix into one df
merge_beta <- 
  merge_beta %>% 
  left_join(melt_env,
            by = c(paste0(level, "1"), paste0(level, "2"))) 

## ---- SGDCs decomposition ----
SGDC.decomp(SD = merge_beta$beta.jtu, 
            GD = merge_beta$Fst, 
            FACTOR = merge_beta[,c("lcdist","distenv")])

