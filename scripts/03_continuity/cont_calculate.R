# This script calculates the "continuity" metric from the respective species 
# and genetic diversity metrics

# load ----
sd_table_global <-
  read_csv("./intermediate/sd_table_global.csv",
           show_col_types = FALSE)

gd_table_global <-
  read_csv("./intermediate/gd_table_global.csv",
           show_col_types = FALSE)

# wrangle ----
cont_global <- list()
for(oc in unique(gd_table_global$ocean)){
  sd <-
    sd_table_global %>% 
    filter(ocean == oc)
  
  gd <-
    gd_table_global %>% 
    filter(ocean == oc) 
  
  cont_global[[oc]] <-
    gd %>% 
    dplyr::select(-ocean) %>% 
    left_join(sd,
              by = "family")
}

# scale and calculate ----
for(oc in names(cont_global)){
  
  # scale
  for(i in which(colnames(cont_global[[oc]]) %in% c("richness_site_mean",
                                                    "beta.JAC",
                                                    "richness_ocean",
                                                    "sd_pd",
                                                    "sd_mpd",
                                                    "sd_vpd",
                                                    "Hs",
                                                    "Fst",
                                                    "Ht",
                                                    "gd_pd",
                                                    "gd_mpd",
                                                    "gd_vpd"))){
    
    cont_global[[oc]][,i] <- 
      scales::rescale(cont_global[[oc]][[i]], to = c(0.01,1))
  }
  
  # calculate
  cont_global[[oc]]$cont_rich_alpha <-
    log(cont_global[[oc]]$richness_site_mean/cont_global[[oc]]$Hs)
  
  cont_global[[oc]]$cont_rich_beta <-
    log(cont_global[[oc]]$beta.JAC/cont_global[[oc]]$Fst)
  
  cont_global[[oc]]$cont_rich_gamma <-
    log(cont_global[[oc]]$richness_ocean/cont_global[[oc]]$Ht)
  
  cont_global[[oc]]$cont_pd <-
    log(cont_global[[oc]]$sd_pd/cont_global[[oc]]$gd_pd)
  
  cont_global[[oc]]$cont_mpd <-
    log(cont_global[[oc]]$sd_mpd/cont_global[[oc]]$gd_mpd)
  
  cont_global[[oc]]$cont_vpd <-
    log(cont_global[[oc]]$sd_vpd/cont_global[[oc]]$gd_vpd)
  
}

# merge ----
cont_global <-
  do.call(rbind.data.frame, cont_global)

# export ----
write_csv(cont_global,
          "./intermediate/cont_global.csv")
