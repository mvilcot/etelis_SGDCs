
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
    metricLAB <- gsub('beta.', 'β', metricSD)
    
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
      left_join(data_sites, by = "site") %>% 
      mutate(site = factor(site, levels = levels(data_sites$site)))
      # left_join(data_stations, by = "station")
      
    # plot
    gg_list[[i]] <-
      ggplot(pcoa_coord, aes(x=A1, y=A2, color=site)) +
      geom_point(size = 3, alpha = 0.7) +
      labs(x=labels[1], y=labels[2]) +
      # geom_text_repel(label = pcoa_coord[[level]], bg.color = "grey70", bg.r = 0.02) +
      scale_color_manual(values = color_perso, labels = LABELS) +
      ggtitle(paste0(metricLAB, ' (', comm, ')')) +
      theme_light()

    i=i+1
  }
  
}

### {FIGURE S5} ####
gg_grob <-
  patchwork::wrap_plots(gg_list, tag_level = "new", guides = "collect") +
  plot_annotation(tag_level = "a", tag_prefix = "(", tag_suffix = ")") 
plot(gg_grob)

ggsave(gg_grob, width = 9, height = 7, dpi = 500,
       filename = paste0("results/2_species_diversity/_S5_PCoA_beta_species_diversity_", comm_delin, ".png"))




