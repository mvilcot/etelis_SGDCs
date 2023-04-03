gd_alpha <- read.csv("results/TEST/gd_table_site_FijiSep.csv")
sd_alpha <- read.csv("results/TEST/sd_table_site_depth_category_sepFiji.csv")


## ---- alpha-SGDCs ----

table_alpha <-
  gd_alpha %>% 
  left_join(sd_alpha, by = level)

level = "site"
metricSD = paste0("richness_", level)
metricGD = "Ho"

gg_list <- list()

for (comm in names(list_communities)){
  table_alpha_comm <- 
    table_alpha %>% 
    filter(community == comm)
  
  sgdc_alpha <- summary(lm(table_alpha_comm[[metricGD]] ~ table_alpha_comm[[metricSD]]))
  
  gg_list[[comm]] <- 
    ggplot(table_alpha_comm, aes(.data[[metricSD]], .data[[metricGD]])) +
    geom_point() +
    xlab(paste0("richness_", level, " (", comm, ")")) +
    annotate('text', 
             x=min(table_alpha_comm[[metricSD]]), y=max(table_alpha_comm[[metricGD]]),
             hjust = 0, vjust = 1,
             label=paste0("LM adj.R² = ", sprintf("%.3f", sgdc_alpha$adj.r.squared), "\np = ", sprintf("%.3f", coef(sgdc_alpha)[2,4])))
  
  # ggsave(gg, device = "png", height = 4, width = 6, units = "in",
  #        filename = paste0("results/04_continuity/SGDCs_alpha_all_", level, "_", metricGD, "_", metricSD, "_", comm, ".png"))
  
}

gg_grob <- arrangeGrob(grobs = gg_list, ncol=3)
plot(gg_grob)
ggsave(gg_grob, width = 15, height = 8, 
       filename = paste0("results/TEST/SGDCs_alpha_all_", level, "_", metricGD, "_", metricSD, "_",  comm_delin, ".png"))

