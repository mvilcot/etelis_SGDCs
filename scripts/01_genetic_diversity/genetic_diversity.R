
library(hierfstat)
library(mmod)
# library("radiator")
# library("adegenet")
# library("HierDpart")
# library("betapart")
# library("diveRsity")

## ---- Set parameters ----
# filters
filters <- "missind_callrate0.70_maf0.05"
level <- "site"

# Read genlight  
genlight <- readRDS(paste0("intermediate/01_genetic_diversity/Genlight_Etelis_coruscans_ordered_", filters, ".RDS"))
genlight

# Set population to site/station
genlight@pop <- genlight@other[["ind.metrics"]][[level]]
genlight

# Remove populations with less than two individuals
print(paste("removing", level, ":", names(which(table(genlight@pop) < 2))))
genlight <- gl.drop.pop(genlight, pop.list = names(which(table(genlight@pop) < 2)), recalc = T, mono.rm = T)
genlight



## ---- Convert to other formats ----
# Genind 
genind <- dartR::gl2gi(genlight) # same than adegenet::df2genind
genind

# Genpop class 
genpop <- adegenet::genind2genpop(genind) 
genpop

# Genepop (Rousset)
# genomic_converter(genlight, output = "genepop")
# genepop <- dartR::gl2genepop(genlight, outfile = 'Intermediate/GenepopdartR_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_allsites.genepop_test2.gen', outpath = ".")
# genepop

# Hierfstat format
# ghierfstat <- hierfstat::genind2hierfstat(genind)



## ---- global genetic diversity ----

# basic stats
BS <- hierfstat::basic.stats(genind) # same result
BS
BSo <- BS$overall # over all samples

# Jost Additive framework 
Hs = BSo[["Hs"]]
Ht = BSo[["Ht"]]
Hst = (Ht-Hs)/(1-Hs)
Dst <- Ht - Hs

# Jost Multiplicative framework Jst=Jt*Js
Js = 1/(1-Hs)
Jt = 1/(1-Ht)
Jst = Jt/Js # or 1/(1-Hst)

# # Gst
# # k=4????????????
# Gst     <- (Ht-Hs)/Ht
# Gst_Max <- (k-1)*(1-Hs)/ (k-1+Hs)
# GstP    <- Gst/Gst_Max
# GstPP   <- k*(Ht-Hs)/ ((k*Ht-Hs)*(1-Hs))

# Gst Hedrick
GstPP_hed <- Gst_Hedrick(genind)

# create global table
table_gd_global <- 
  cbind(Ho=BSo["Ho"], Hs, Ht, Hst, Dst, 
        Js, Jt, Jst,
        Fst=BSo["Fst"], Fis=BSo["Fis"],
        Dstp=BSo["Dstp"], Dest=BSo["Dest"],
        # Gst, GstP, GstPP,
        Gstpp_hed=GstPP_hed$global) %>% 
  as_tibble()


## ---- alpha gd by location ----
table_gd_alpha <-
  data.frame(Hs = colMeans(BS$Hs, na.rm = T))

# table_gd_alpha_site2 <- 
#   as.data.frame(BS$Hs) %>% 
#   rownames_to_column("loci") %>%
#   pivot_longer(-loci, names_to = "site") %>% 
#   group_by(site) %>% 
#   summarise(Hs=mean(value, na.rm=T), .groups = "keep")

# table_gd_alpha_site3 <- 
#   data.frame(Hs = apply(BS$Hs, MARGIN = 2, FUN = mean, na.rm = T)) # mean by population


## ---- beta gd pairwise ----

# Fst
Fst_pair <- genet.dist(genind, method = "WC84")

# G"st
# Gst_Hedrick(genind) # Gst global = 0.045
GstPP_hed_pair <- pairwise_Gst_Hedrick(genind)


list_gd_beta_pair <-
  list(Fst_pair, GstPP_hed_pair)




## ---- export ----
write.csv(table_gd_global, "results/02_species_diversity/sd_table_global", level, ".csv", row.names = FALSE)
write.csv(table_gd_alpha, "results/02_species_diversity/sd_table_", level, ".csv", row.names = FALSE)
saveRDS(list_gd_beta_site, "results/02_species_diversity/sd_list_pairwise_", level, ".RDS")



## ---- DRAFTS ----

# ## TEST PAIRWISE (diveRsity)
# test <- diffCalc(infile = "Intermediate/GenepopRadiator_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_sites_noCocos.gen", 
#          outfile = "test", fst = T, pairwise = T)
# 
# plot(test$pairwise$D ~ test$pairwise$GGst)
# 
# 
# test2 <- fastDivPart(infile = "Intermediate/GenepopRadiator_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_sites_noCocos.gen", outfile = "myresults", gp = 3, 
#                  pairwise = T, plot = TRUE)




# ## Jaccard
# GDbetaJAC <- HierJd("Intermediate/test.gen", ncode = 3, r = 1, nreg = 1)
# saveRDS(GDbetaJAC, "Results/test_HierJd_betaGD_stations_noinf2.RDS")
# # !!!!!!!!!!!!!!!!! Error in HierJd("Intermediate/test.gen", ncode = 3) : !!!!!!!!!!!! ####
# # argument "r" is missing, with no default
# # In addition: There were 27 warnings (use warnings() to see them)
# 
# ########### ICI ################
# ## Jaccard Betapart
# PAalleles <- as.data.frame(genind@tab)
# PAallelesTEST <- PAalleles
# PAallelesTEST[PAallelesTEST > 0] <- 1 
# 
# PAallelesTEST$site <- genind@pop
# PAallelesTESTagg <- aggregate(PAallelesTEST[,-ncol(PAallelesTEST)], by=list(Site=PAallelesTEST$site), FUN=max) ##by site
# 
# 
# test <- beta.multi(PAallelesTEST[,-ncol(PAallelesTEST)], index.family="jaccard")
# 




