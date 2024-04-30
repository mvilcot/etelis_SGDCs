
# ---- read SNPs dataset ----

# parameters
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "allsites"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)

genlight 


# ---- convert formats ----

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

gmatrixPAloc <- as.data.frame(gmatrixPA)
gmatrixPAloc[[level]] <- genlight@pop
gmatrixPAloc <- 
  aggregate(gmatrixPAloc[,-ncol(gmatrixPAloc)], by=list(location=gmatrixPAloc[[level]]), FUN=max) %>% ##by site
  column_to_rownames("location")



# ---- mean genetic diversity ----
# basic stats
BS <- basic.stats(genind)
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
GstPP.hed <- mmod::Gst_Hedrick(genind)

# Jost D (Jost 2008)
JostD <- mmod::D_Jost(genind)

# Fst population-specific (Weir and Goudet 2017)
popFst <- hierfstat::betas(ghierfstat, nboot=1000)

# Jaccard
jac_multi <- betapart::beta.multi(gmatrixPAloc, index.family="jaccard")

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
  dplyr::as_tibble()



# ---- alpha gd by location ----
gd_alpha <-
  data.frame(Hs = colMeans(BS$Hs, na.rm = T),
             Ho = colMeans(BS$Ho, na.rm = T)) %>% 
  tibble::rownames_to_column(level) %>% 
  dplyr::left_join(
    data.frame(popFst.WG = popFst$betaiovl) %>% 
      tibble::rownames_to_column(level), 
    by = level)



# ---- beta gd pairwise ----
# Fst
Fst_pair <- hierfstat::genet.dist(genind, method = "WC84")

# Hedrick G"st
GstPP.hed_pair <- mmod::pairwise_Gst_Hedrick(genind)

# Jost D
JostD_pair <- mmod::pairwise_D(genind)

# Jaccard
jac_pair <- betapart::beta.pair(gmatrixPAloc, index.family="jaccard")

# put into list
list_gd_beta_pair <-
  list(Fst = Fst_pair, 
       GstPP.hed = GstPP.hed_pair,
       D.Jost = JostD_pair,
       jtu = jac_pair$beta.jtu,
       jac = jac_pair$beta.jac,
       jne = jac_pair$beta.jne)



# ---- export ----
BS %>% saveRDS(paste0("intermediate/1_genetic_diversity/basic_stats_", level, ".RDS"))
gd_global %>% write.csv(paste0("results/1_genetic_diversity/gd_table_global_", level, ".csv"), row.names = F, quote = F)
gd_alpha %>% write.csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"), row.names = F, quote = F)
list_gd_beta_pair %>% saveRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))


