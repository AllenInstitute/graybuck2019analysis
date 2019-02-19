library(lowcat)
library(scrattch.io)
library(scrattch.vis)
library(feather)
library(dplyr)

# River plot
source("sankey_functions.R")

anno <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/rnaseq_analysis/doc/paper/20180622_mscRE_mapping/mscre4_visp_mapping.csv")

anno <- anno %>%
  rename_("donor" = "donor_label") %>%
  annotate_cat("donor")

nodes <- make_plot_nodes(make_group_nodes(anno, c("donor","pred_cluster")), pad = 0.2)

right_nodes <- nodes %>%
  filter(group == "pred_cluster") %>%
  mutate(name = paste0(name, " (",n,")"))

river_plot <- build_river_plot(anno,
                 c("donor","pred_cluster"),
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

ggsave("mscRE4_scRNAseq_mapping_river.pdf",
       width = 3, height = 3,
       useDingbats = F)
