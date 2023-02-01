
library("ggplot2")
library("vegan")
library("ecodist")


# ----------- Setup data -------------------------------------------------------

sitesmeta <- read.csv("Data/sites_coordinates.csv")

## Read GD and SD distance matrix
SDbeta <- readRDS("Results/SD_beta_stations_betapair.RDS")
SDbeta <- as.matrix(SDbeta$beta.jtu)
GDbeta <- readRDS("Results/GD_beta_stations_noinf2_pairwise_GstppHedrick.RDS")
GDbeta <- as.matrix(GDbeta)

#### !!!! CHANGE SD row names #######
rownames(SDbeta) <- colnames(SDbeta) <- 
  sitesmeta[match(rownames(SDbeta), sitesmeta$Station), ]$Full.name


## Keep only sampling sites present in both GD and SD
SDbeta <- SDbeta[rownames(SDbeta) %in% rownames(GDbeta), 
                 colnames(SDbeta) %in% colnames(GDbeta)]
GDbeta <- GDbeta[rownames(GDbeta) %in% rownames(SDbeta), 
                 colnames(GDbeta) %in% colnames(SDbeta)]

## Order rows alphabetically
SDbeta <- SDbeta[order(rownames(SDbeta)), order(colnames(SDbeta))]
GDbeta <- GDbeta[order(rownames(GDbeta)), order(colnames(GDbeta))]


## Merge two distance matrix into one data.frame
meltSD <- t(combn(colnames(SDbeta), 2))
meltSD <- data.frame(meltSD, jtu=SDbeta[meltSD])

meltGD <- t(combn(colnames(GDbeta), 2))
meltGD <- data.frame(meltGD, Gstpp=GDbeta[meltGD])

beta <- merge(meltSD, meltGD, all = T)
beta$sites <- paste0(beta$X1, "-", beta$X2)



# ----------- beta-SGDCs -------------------------------------------------------

## Mantel test
vegan::mantel(as.dist(SDbeta), as.dist(GDbeta))

## MRM
MRM(formula = jtu ~ Gstpp, data = beta, nperm = 999)

## Plot
ggplot(data = beta, aes(jtu, Gstpp, label=sites)) +
  geom_point()
ggsave("Results/07_SGDCs/Etelis_SGDCs.png")



# ----------- beta-SGDCs without Seychelles ------------------------------------

patt = "Seychelles"
betaSUB <- beta[!grepl(patt, beta$sites),]

## Mantel test
SDbetaSUB <- SDbeta[!grepl(patt, rownames(SDbeta)),
                    !grepl(patt, colnames(SDbeta))]
GDbetaSUB <- GDbeta[!grepl(patt, rownames(GDbeta)),
                    !grepl(patt, colnames(GDbeta))]
vegan::mantel(as.dist(SDbetaSUB), as.dist(GDbetaSUB))

## MRM
MRM(formula = jtu ~ Gstpp, data = betaSUB, nperm = 999)

## Plot
ggplot(data = betaSUB, aes(jtu, Gstpp, label=sites)) +
  geom_point() 
ggsave("Results/07_SGDCs/Etelis_SGDCs_noSeychelles.png")




# ----------- beta-SGDCs without Hawaii and Seychelles ------------------------------------

patt = "Hawaii|Seychelles"
betaSUB <- beta[!grepl(patt, beta$sites),]

## Mantel test
SDbetaSUB <- SDbeta[!grepl(patt, rownames(SDbeta)),
                    !grepl(patt, colnames(SDbeta))]
GDbetaSUB <- GDbeta[!grepl(patt, rownames(GDbeta)),
                    !grepl(patt, colnames(GDbeta))]
vegan::mantel(as.dist(SDbetaSUB), as.dist(GDbetaSUB))

## MRM
MRM(formula = jtu ~ Gstpp, data = betaSUB, nperm = 999)

## Plot
ggplot(data = betaSUB, aes(jtu, Gstpp, label=sites)) +
  geom_point() 
ggsave("Results/07_SGDCs/Etelis_SGDCs_noSeychellesHawaii.png")






# ----------- beta-SGDCs without Hawaii and Seychelles ------------------------------------

patt = "Hawaii"
betaSUB <- beta[!grepl(patt, beta$sites),]

## Mantel test
SDbetaSUB <- SDbeta[grepl(patt, rownames(SDbeta)),
                    grepl(patt, colnames(SDbeta))]
GDbetaSUB <- GDbeta[grepl(patt, rownames(GDbeta)),
                    grepl(patt, colnames(GDbeta))]
vegan::mantel(as.dist(SDbetaSUB), as.dist(GDbetaSUB))

## MRM
MRM(formula = jtu ~ Gstpp, data = betaSUB, nperm = 999)

## Plot
ggplot(data = betaSUB, aes(jtu, Gstpp, label=sites)) +
  geom_point() 
ggsave("Results/07_SGDCs/Etelis_SGDCs_noSeychellesHawaii.png")

