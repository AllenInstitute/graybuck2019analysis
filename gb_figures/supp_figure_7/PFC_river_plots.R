library(scrattch.vis)
library(scrattch.io)
library(ggplot2)
library(dplyr)
library(purrr)
options(stringsAsFactors = F)

# Cusanovich 2018 and Graybuck 2019
mapping <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2019-01-11_Graybuck_and_GSE111586/2e4_tss_phenograph_scatac_visp_marker_correlation.csv",
                    row.names = 1) %>%
  select(-x,-y)

cu_meta <- read.table("cell_metadata.tissue_freq_filtered.txt", sep = "\t",
                      header = TRUE) %>%
  mutate(barcode = substr(cell, 1, 32))
names(cu_meta)[names(cu_meta) == "id"] <- "cu_type"

cu_samples <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2019-01-11_GSE111586_analysis/all_samples.csv")

cu_type_anno <- cu_meta %>%
  select(cu_type) %>%
  unique() %>%
  annotate_cat(cu_type)

anno <- mapping %>%
  left_join(cu_samples) %>%
  left_join(cu_meta) %>%
  filter(!is.na(barcode)) %>%
  mutate(cu_type_label = cu_type) %>%
  left_join(cu_type_anno)

# Filter cu_types that don't have > 80% matching a single class
anno <- anno %>%
  mutate(class_id = ifelse(pred_cluster_id %in% 1:60,
                           1, ifelse(pred_cluster_id < 117, 2, 3)))

anno <- anno %>%
  group_by(cu_type_id) %>%
  mutate(keep_id = max(table(class_id)) > 0.8*n()) %>%
  ungroup()

anno <- anno %>%
  filter(keep_id)

# Update the cu_type_id based on mapping to pred_cluster_id
new_order <- anno %>%
  group_by(cu_type_id) %>%
  summarise(mean_cl = mean(pred_cluster_id)) %>%
  arrange(mean_cl) %>%
  mutate(cu_type_id2 = 1:n())

anno <- anno %>%
  left_join(new_order) %>%
  mutate(cu_type_id = cu_type_id2)

river <- build_river_plot(anno,
                          c("cu_type","pred_cluster"),
                          min_link_size = 0.05)

ggsave("cu_to_tasic_river.pdf",
       river,
       width = 2,
       height = 4.5)
