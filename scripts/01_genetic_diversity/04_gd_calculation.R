
# library(HierDpart)


## ---- set parameters ----
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

# exclusion_list_site <- c("Cocos_Keeling")
# exclusion_list_station <- c("Cocos (Keeling)", 
#                             "New Caledonia (Antigonia)",
#                             "New Caledonia (Jumeau)", 
#                             "New Caledonia (Poum)")




## ---- convert to other formats ----
# Genind 
genind <- dartR::gl2gi(genlight) # same than adegenet::df2genind
genind

# # Genpop class
# genpop <- adegenet::genind2genpop(genind)
# genpop
# 
# # Genepop (Rousset)
# genomic_converter(genlight, output = "genepop")
# genepop <- dartR::gl2genepop(genlight, outfile = paste0('Intermediate/Genepop_Etelis_coruscans_ordered_', filters, '.genepop.gen', outpath = "."))
# genepop
# 
# # Hierfstat format
# ghierfstat <- hierfstat::genind2hierfstat(genind)
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


## ---- global genetic diversity ----

# basic stats
BS <- hierfstat::basic.stats(genind) # same result
BS
BSo <- BS$overall # over all samples

# Jost Additive framework Dst = Ht - Hs
Hs = BSo[["Hs"]]
Ht = BSo[["Ht"]]
Hst = (Ht-Hs)/(1-Hs)
Dst <- Ht - Hs

# Jost Multiplicative framework Jst = Jt / Js
Js = 1/(1-Hs)
Jt = 1/(1-Ht)
Jst = Jt/Js # or 1/(1-Hst)

# Hedrick Gst"
GstPP.hed <- Gst_Hedrick(genind)

# Jost D
JostD <- D_Jost(genind)

# Confidence interval 
bs <- chao_bootstrap(genind) 
bs_D <- summarise_bootstrap(bs, JostD) 
bias <- bs.D$summary.global.het[1]- obs.D$global.het 
bs.D$summary.global.het- bias


# create global table
gd_global <-
  cbind(Ho = BSo["Ho"], Hs, Ht, Hst, Dst,
        Js, Jt, Jst,
        Fst = BSo["Fst"], 
        Fis = BSo["Fis"],
        Dstp = BSo["Dstp"], 
        Dest = BSo["Dest"],
        Gstpp.hed = GstPP.hed$global,
        JostD = JostD$global.het) %>%
  as_tibble()




## ---- alpha gd by location ----
gd_alpha <-
  data.frame(Hs = colMeans(BS$Hs, na.rm = T),
             Ho = colMeans(BS$Ho, na.rm = T)) %>% 
  rownames_to_column(var = level)



## ---- beta gd pairwise ----

# Fst
Fst_pair <- genet.dist(genind, method = "WC84")

# G"st
GstPP.hed_pair <- pairwise_Gst_Hedrick(genind)

# Jost D
JostD_pair <- pairwise_D(genind)


# put into list
list_gd_beta_pair <-
  list(Fst = Fst_pair, 
       GstPP.hed = GstPP.hed_pair,
       JostD = JostD_pair)




## ---- export ----
saveRDS(BS, "intermediate/01_genetic_diversity/basic_stats.RDS")
write.csv(gd_global, paste0("results/01_genetic_diversity/gd_table_global_", level, ".csv"), row.names = F, quote = F)
write.csv(gd_alpha, paste0("results/01_genetic_diversity/gd_table_", level, ".csv"), row.names = F, quote = F)
saveRDS(list_gd_beta_pair, paste0("results/01_genetic_diversity/gd_list_pairwise_", level, ".RDS"))





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
# 
# 
# 
# 
# ## Jaccard
# GDbetaJAC <- HierJd("Intermediate/test.gen", ncode = 3, r = 1, nreg = 1)
# saveRDS(GDbetaJAC, "Results/test_HierJd_betaGD_stations_noinf2.RDS")
# # !!!!!!!!!!!!!!!!! Error in HierJd("Intermediate/test.gen", ncode = 3) : !!!!!!!!!!!!
# # argument "r" is missing, with no default
# # In addition: There were 27 warnings (use warnings() to see them)
# 
# 
# 
# ## TEST Jaccard Betapart
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
# 
# # Jaccard
# JacHier <- HierJd("_archive/-36_radiator_genomic_converter_20221115@1129/radiator_data_20221115@1129_genepop.gen", 
#                   ncode=3, nreg=1, r=1)
# print(JacHier)
# 
# 
# library(jacpop)
# jac <- generate_pw_jaccard(geno = test)




# gd_alpha_site2 <-
#   as.data.frame(BS$Hs) %>%
#   rownames_to_column("loci") %>%
#   pivot_longer(-loci, names_to = "site") %>%
#   group_by(site) %>%
#   summarise(Hs=mean(value, na.rm=T), .groups = "keep")

# gd_alpha_site3 <-
#   data.frame(Hs = apply(BS$Hs, MARGIN = 2, FUN = mean, na.rm = T)) # mean by population
