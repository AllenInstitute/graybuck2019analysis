library(lowcat)
library(scrattch.io)
library(scrattch.vis)
library(feather)
library(dplyr)
options(stringsAsFactors = FALSE)

# River plot
source("sankey_functions.R")

anno <- read.csv("mscre_visp_mapping_for_plots.csv")

anno %>%
  group_by(mscre_label) %>%
  summarise(n_cells = n(),
            n_lt_0.8 = sum(cor <= 0.8))

anno <- anno %>%
  filter(cor > 0.8)

mscres <- c("mscRE4-FlpO","mscRE10-FlpO","mscRE16-FlpO")

walk(mscres, function(x) {
  anno <- anno %>% filter(mscre_label == x)
  nodes <- make_plot_nodes(make_group_nodes(anno, c("mscre","pred_cluster")), pad = 0.2)
  
  right_nodes <- nodes %>%
    filter(group == "pred_cluster") %>%
    mutate(name = paste0(name, " (",n,")"))
  
  river_plot <- build_river_plot(anno,
                                 c("mscre","pred_cluster"),
                                 pad = 0.2,
                                 fill_group = "pred_cluster") +
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
  nodes <- make_plot_nodes(make_group_nodes(anno, c("mscre","pred_subclass","pred_cluster")), pad = 0.2)
  
  right_nodes <- nodes %>%
    filter(group == "pred_cluster") %>%
    mutate(name = paste0(name, " (",n,")"))
  
  river_plot <- build_river_plot(anno,
                                 c("mscre","pred_subclass","pred_cluster"),
                                 pad = 0.2,
                                 fill_group = "pred_subclass") +
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

