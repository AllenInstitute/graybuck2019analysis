library(lowcat)
library(scrattch.io)
library(scrattch.vis)
library(feather)
library(dplyr)
options(stringsAsFactors = FALSE)

# River plot
source("sankey_functions.R")

anno <- read.csv("rnaseq_sample_metadata_with_mapping.csv")

anno %>%
  group_by(experiment) %>%
  summarise(n_cells = n(),
            n_lt_0.8 = sum(prob <= 0.8))

anno <- anno %>%
  filter(prob > 0.8) %>%
  annotate_cat(experiment) %>%
  annotate_cat(external_donor_name)

mscres <- c("mscRE16-EGFP_st_left",
            "mscRE16-EGFP_st_right",
            "mscRE4-EGFP_st_left",
            "mscRE4-EGFP_st_right")

walk(mscres, function(x) {
  anno <- anno %>% filter(experiment_label == x)
  nodes <- make_plot_nodes(make_group_nodes(anno, c("external_donor_name","cluster")), pad = 0.2)
  
  right_nodes <- nodes %>%
    filter(group == "cluster") %>%
    mutate(name = paste0(name, " (",n,")"))
  
  river_plot <- build_river_plot(anno,
                                 c("external_donor_name","cluster"),
                                 pad = 0.2,
                                 fill_group = "cluster") +
    geom_text(data = right_nodes,
              aes(x = xmax,
                  y = (ymin + ymax)/2,
                  label = name,
                  color = color),
              hjust = 0,
              size = 16/6) +
    scale_color_identity()
  
  ggsave(paste0(x,"_scRNAseq_mapping_river.pdf"),
         river_plot,
         width = 3, height = 3,
         useDingbats = F)
})


walk(mscres, function(x) {
  anno <- anno %>% filter(mscre_label == x)
  nodes <- make_plot_nodes(make_group_nodes(anno, c("mscre","subclass")), pad = 0.2)
  
  right_nodes <- nodes %>%
    filter(group == "subclass") %>%
    mutate(name = paste0(name, " (",n,")"))
  
  river_plot <- build_river_plot(anno,
                                 c("mscre","subclass"),
                                 pad = 0.2,
                                 fill_group = "subclass") +
    geom_text(data = right_nodes,
              aes(x = xmax,
                  y = (ymin + ymax)/2,
                  label = name,
                  color = color),
              hjust = 0,
              size = 16/6) +
    scale_color_identity()
  
  ggsave(paste0(x,"_scRNAseq_subclass_mapping_river.pdf"),
         river_plot,
         width = 3, height = 3,
         useDingbats = F)
})

