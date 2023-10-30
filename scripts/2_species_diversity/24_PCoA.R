
# ---- read distance matrix ----
level = "site"
comm_delin = "taxonomy"

sd_mat <- readRDS(paste0("results/2_species_diversity/sd_list_pairwise_", level, "_", comm_delin, ".RDS"))






# ---- PCoA ----

gg_list <- list()
i=1
for (metricSD in c("beta.jac", "beta.jtu")){
  # for (comm in names(sd_mat)){
  for (comm in names(sd_mat)[c(2,4)]){
      # read distance matrix
    mat_SDbeta <- sd_mat[[comm]][[metricSD]]
    
    # run pcoa
    pcoa <- dudi.pco(quasieuclid(mat_SDbeta), scannf=FALSE, nf=8)

    # eigenvalues
    percent_explained <- pcoa$eig / sum(pcoa$eig) * 100
    
    # set labels
    pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
    labels <- c(glue("PCo Axis 1 ({pretty_pe[1]}%)"),
                glue("PCo Axis 2 ({pretty_pe[2]}%)"))
    
    pcoa_coord <-
      pcoa$li %>%
      rownames_to_column(level) %>% 
      left_join(data_sites, by = "site")
      # left_join(data_stations, by = "station")
      
    # plot
    # gg_list[[i]] <-
      ggplot(pcoa_coord, aes(x=A1, y=A2, color=site)) +
      geom_point() +
      labs(x=labels[1], y=labels[2]) +
      geom_text_repel(label = pcoa_coord[[level]], bg.color = "grey70", bg.r = 0.02) +
      scale_color_manual(values = color_perso) +
      ggtitle(comm, subtitle = metricSD) +
      theme_light() +
      theme(legend.position="none")
    
    i=i+1
  }
  
}

gg_grob <-
  patchwork::wrap_plots(gg_list, tag_level = "new") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")") 
plot(gg_grob)

ggsave(gg_grob, width = 10, height = 10, 
       filename = paste0("results/2_species_diversity/PCoA_beta_species_diversity_", comm_delin, "-2.png"))



