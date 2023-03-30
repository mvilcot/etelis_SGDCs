
## ---- parameters ----
# filters
filters = "missind_callrate0.70_maf0.05"
level = "site"

# read genlight  
genlight <- 
  read.genlight(filters, level, removeless2ind = TRUE)





## ---- convert to other formats ----

# genind 
genind <- dartR::gl2gi(genlight)

# hierfstat format
ghierfstat <- hierfstat::genind2hierfstat(genind)

# SNP presence/absence lfmm (package LEA, SilicoDArT)
geno <- dartR::gl2geno(genlight)

# matrix
gmatrix <- as.matrix(genlight) # one column by SNP
gmatrixPA <- gmatrix[ , colSums(is.na(gmatrix))==0]
gmatrixPA[gmatrixPA >= 1] <- 1
# gmatrix2 <- as.matrix(genind@tab) # one column by allele
# gmatrix2PA <- gmatrix2[ , colSums(is.na(gmatrix2))==0]
# gmatrix2PA[gmatrix2PA >= 1] <- 1

gmatrixPAloc <- as.data.frame(gmatrixPA)
gmatrixPAloc[[level]] <- genlight@pop
gmatrixPAloc <- 
  aggregate(gmatrixPAloc[,-ncol(gmatrixPAloc)], by=list(location=gmatrixPAloc[[level]]), FUN=max) %>% ##by site
  column_to_rownames("location")



## ---- tess3r ----
# https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html

library(tess3r)





## ---- mean genetic diversity ----

# basic stats
BS <- hierfstat::basic.stats(genind)
BSo <- BS$overall

# Jost additive framework Dst = Ht - Hs
Hs = BSo[["Hs"]]
Ht = BSo[["Ht"]]
Hst = (Ht-Hs)/(1-Hs)
Dst <- Ht - Hs

# Jost multiplicative framework Jst = Jt / Js
Js = 1/(1-Hs)
Jt = 1/(1-Ht)
Jst = Jt/Js # or 1/(1-Hst)

# Hedrick Gst" (Meirmans and Hedrick 2011)
GstPP.hed <- Gst_Hedrick(genind)

# Jost D (Jost 2008)
JostD <- D_Jost(genind)

# Fst population-specific (Weir and Goudet 2017)
popFst <- betas(ghierfstat, nboot=1000)

# # Confidence interval
# bs <- chao_bootstrap(genind)
# bs_D <- summarise_bootstrap(bs, JostD)
# bias <- bs.D$summary.global.het[1]- obs.D$global.het
# bs.D$summary.global.het- bias

# Jaccard
### either pairwise by individual, then average
# jac <- beta.multi(gmatrixPA, index.family="jaccard")
# jac2 <- beta.pair(gmatrix2PA, index.family="jaccard")
# 
# jacSIM <- generate_pw_jaccard(geno = gmatrix, 
#                            pop.label = genlight$pop,
#                            plot_it = FALSE)

### or PA by site first, then pairwise jaccard
jac_multi <- beta.multi(gmatrixPAloc, index.family="jaccard")

# create global table
gd_global <-
  cbind(Ho = BSo["Ho"], Hs, Ht, Hst, Dst,
        Js, Jt, Jst,
        Fst = BSo["Fst"], 
        Fis = BSo["Fis"],
        Dstp = BSo["Dstp"], 
        Dest = BSo["Dest"],
        Gstpp.hed = GstPP.hed$global,
        D.Jost = JostD$global.het,
        popFst.WG = popFst$betaW,
        jtu = jac_multi$beta.JTU,
        jac = jac_multi$beta.JAC,
        jne = jac_multi$beta.JNE) %>%
  as_tibble()




## ---- alpha gd by location ----
gd_alpha <-
  data.frame(Hs = colMeans(BS$Hs, na.rm = T),
             Ho = colMeans(BS$Ho, na.rm = T)) %>% 
  rownames_to_column(level) %>% 
  left_join(
    data.frame(popFst.WG = popFst$betaiovl) %>% 
              rownames_to_column(level), 
            by = level)




## ---- beta gd pairwise ----

# Fst
Fst_pair <- genet.dist(genind, method = "WC84")

# Hedrick G"st
GstPP.hed_pair <- pairwise_Gst_Hedrick(genind)

# Jost D
JostD_pair <- pairwise_D(genind)

# Jaccard
### either PA by site first, then pairwise jaccard
jac_pair <- beta.pair(gmatrixPAloc, index.family="jaccard")

### or pairwise by individual, then average
jac_pair <- beta.pair(gmatrixPA, index.family="jaccard")



# put into list
list_gd_beta_pair <-
  list(Fst = Fst_pair, 
       GstPP.hed = GstPP.hed_pair,
       D.Jost = JostD_pair,
       jtu = jac_pair$beta.jtu,
       jac = jac_pair$beta.jac,
       jne = jac_pair$beta.jne)


## ---- export ----
saveRDS(BS, paste0("intermediate/01_genetic_diversity/basic_stats_", level, ".RDS"))
write.csv(gd_global, paste0("results/01_genetic_diversity/gd_table_global_", level, ".csv"), row.names = F, quote = F)
write.csv(gd_alpha, paste0("results/01_genetic_diversity/gd_table_", level, ".csv"), row.names = F, quote = F)
saveRDS(list_gd_beta_pair, paste0("results/01_genetic_diversity/gd_list_pairwise_", level, ".RDS"))



## ---- DRAFTS ----

## ---- Pairwise (diveRsity) ---- #
# test <- diffCalc(infile = "Intermediate/GenepopRadiator_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_sites_noCocos.gen", 
#          outfile = "test", fst = T, pairwise = T)
# 
# plot(test$pairwise$D ~ test$pairwise$GGst)
# 
# 
# test2 <- fastDivPart(infile = "Intermediate/GenepopRadiator_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_sites_noCocos.gen", outfile = "myresults", gp = 3, 
#                  pairwise = T, plot = TRUE)



## ---- Jaccard ---- #
# ## HierJd 
# GDbetaJAC <- HierJd("Intermediate/test.gen", ncode = 3, r = 1, nreg = 1)
# saveRDS(GDbetaJAC, "Results/test_HierJd_betaGD_stations_noinf2.RDS")
# # !!!!!!!!!!!!!!!!! Error in HierJd("Intermediate/test.gen", ncode = 3) : !!!!!!!!!!!!
# # argument "r" is missing, with no default
# # In addition: There were 27 warnings (use warnings() to see them)
# 
# 
# 
## TEST Jaccard Betapart
PAalleles <- as.matrix(as.data.frame(genind@tab))
PAallelesTEST <- PAalleles[ , colSums(is.na(PAalleles))==0]
PAallelesTEST[PAallelesTEST > 0] <- 1
test <- beta.pair(PAallelesTEST, index.family="jaccard")

PAallelesTEST$site <- genind@pop
PAallelesTESTagg <- aggregate(PAallelesTEST[,-ncol(PAallelesTEST)], by=list(Site=PAallelesTEST$site), FUN=max) ##by site
test <- beta.multi(PAallelesTEST[,-ncol(PAallelesTEST)], index.family="jaccard")
# 
# 
# Jaccard
# library(HierDpart)
# JacHier <- HierJd("_archive/-36_radiator_genomic_converter_20221115@1129/radiator_data_20221115@1129_genepop.gen",
#                   ncode=3, nreg=1, r=1)
# print(JacHier)
# 




## ---- Compute mean alpha by site ---- #
# gd_alpha_site2 <-
#   as.data.frame(BS$Hs) %>%
#   rownames_to_column("loci") %>%
#   pivot_longer(-loci, names_to = "site") %>%
#   group_by(site) %>%
#   summarise(Hs=mean(value, na.rm=T), .groups = "keep")
# 
# gd_alpha_site3 <-
#   data.frame(Hs = apply(BS$Hs, MARGIN = 2, FUN = mean, na.rm = T)) # mean by population
# 
# 
# 
# 
## ---- convert to other formats ---- #
# # Genpop class
# genpop <- adegenet::genind2genpop(genind)
# genpop
# 
# # Genepop (Rousset)
# genomic_converter(genlight, output = "genepop")
# genepop <- dartR::gl2genepop(genlight, outfile = paste0('Intermediate/Genepop_Etelis_coruscans_ordered_', filters, '.genepop.gen', outpath = "."))
# genepop
# 
# # data frame format
# gdf <- genind2df(genind)
# gdf2 <- 
#   gdf %>% 
#   select(-pop)
# 
# alleles <- 
#   read.table(text = gsub("/", "-", colnames(gdf2)), sep = "-", as.is = TRUE) %>%# get base variant for each SNP
#   select(3, 4) %>%
#   cbind(colnames(gdf2), .) %>%
#   dplyr::rename(name = 1, var1 = 2, var2 = 3)
#   
# for (i in 1:nrow(alleles)){
#   allele1 <- alleles[i, "var1"]
#   allele2 <- alleles[i, "var2"]
#   
#   gdf2[!is.na(gdf[,i]) & gdf2[,i] == paste0(allele1, allele1), i] <- 0
#   gdf2[!is.na(gdf[,i]) & gdf2[,i] == paste0(allele1, allele2), i] <- 1
#   gdf2[!is.na(gdf[,i]) & gdf2[,i] == paste0(allele2, allele1), i] <- 1
#   gdf2[!is.na(gdf[,i]) & gdf2[,i] == paste0(allele2, allele2), i] <- 2
#   
# }
# test <- as.data.frame(sapply(gdf2, as.numeric))
# rownames(test) <- rownames(gdf2)
# 
# 
# 
# # plink
# genlight$chromosome <- as.factor("1")
# gl2plink(genlight, bed_file = F, 
#          outpath = getwd(), outfile = "test.plink")


