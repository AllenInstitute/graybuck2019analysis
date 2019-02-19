library(ggplot2)
library(scrattch.vis)
library(scrattch.io)
options(stringsAsFactors = F)

# Cusanovich 2018 and Graybuck 2019
load("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2019-01-11_Graybuck_and_GSE111586/f1e4_e25_c10_ds1e4_x1e3_jaccard.rda")

mapping <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2019-01-11_Graybuck_and_GSE111586/2e4_tss_phenograph_scatac_visp_marker_correlation.csv",
                    row.names = 1) %>%
  select(-x,-y)

cu_meta <- read.table("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE111586/online_metadata/cell_metadata.tissue_freq_filtered.txt", sep = "\t",
                      header = TRUE) %>%
  mutate(barcode = substr(cell, 1, 32))

names(cu_meta)[names(cu_meta) == "id"] <- "cu_type"

cell_label_anno <- read.csv("cell_label_colors.csv")

hill_coords <- read.table("andrewhill_figure5g_coordinates.txt", header = TRUE) %>%
  mutate(barcode = substr(cell, 1, 32))

names(hill_coords) <- c("cell","tsne1_ah","tsne2_ah","barcode")

cu_samples <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2019-01-11_GSE111586_analysis/all_samples.csv")
cu_samples <- cu_samples %>%
  mutate(dataset_full = "Cusanovich_2018",
         dataset_short = "cu") %>%
  left_join(cu_meta, by = "barcode") %>%
  left_join(hill_coords, by = "barcode")

plot_df <- tsne_df %>%
  left_join(cu_samples) %>%
  mutate(cu_type_label = ifelse(is.na(cu_type), "Graybuck", cu_type)) %>%
  left_join(mapping)

cu_type_anno <- cu_meta %>%
  select(cu_type) %>%
  unique() %>%
  annotate_cat(cu_type)

cu_type_anno <- rbind(cu_type_anno, 
                      data.frame(cu_type_label = "Graybuck",
                                 cu_type_id = 0,
                                 cu_type_color = "#808080"))

plot_df <- plot_df %>%
  left_join(cu_type_anno)

cu_meta <- cu_meta %>%
  filter(tissue == "PreFrontalCortex") %>%
  left_join(hill_coords) %>%
  mutate(cu_type_label = cu_type) %>%
  left_join(cu_type_anno) %>%
  mutate(downstream = ifelse(barcode %in% plot_df$barcode,TRUE, FALSE))

cu_tsne <- ggplot() +
  geom_point(data = cu_meta %>% filter(downstream == FALSE),
             aes(x = tsne1_ah,
                 y = tsne2_ah,
                 color = cu_type_color),
             size = 0.1,
             alpha = 0.15) +
  geom_point(data = cu_meta %>% filter(downstream == TRUE),
             aes(x = tsne1_ah,
                 y = tsne2_ah,
                 color = cu_type_color),
             size = 0.1) +
  scale_color_identity() +
  theme_classic(base_size = 5)

ggsave("cu_tsne.pdf",
       cu_tsne,
       width = 2, height = 2,
       useDingbats = FALSE)

joint_tsne <- ggplot() +
  geom_point(data = plot_df,
             aes(x = x,
                 y = y,
                 color = cu_type_color),
             size = 0.1) +
  scale_color_identity() +
  theme_classic(base_size = 5)

ggsave("joint_tsne_cu_colors.pdf",
       joint_tsne,
       width = 2, height = 2,
       useDingbats = FALSE)

joint_tsne_cluster <- ggplot() +
  geom_point(data = plot_df,
             aes(x = x,
                 y = y,
                 color = pred_cluster_color),
             size = 0.1) +
  scale_color_identity() +
  theme_classic(base_size = 5)

ggsave("joint_tsne_cluster_colors.pdf",
       joint_tsne_cluster,
       width = 2, height = 2,
       useDingbats = FALSE)


# cu-only lowcat tsne
cu_only_mapping <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2019-01-11_GSE111586_analysis/2e4_tss_phenograph_scatac_visp_marker_correlation.csv",
                            row.names = 1)

sample_to_cu <- plot_df %>%
  select(sample_id, cu_type_label, cu_type_id, cu_type_color)

cu_only_df <- left_join(cu_only_mapping,
                        sample_to_cu)

cu_only_tsne <- ggplot() +
  geom_point(data = cu_only_df,
             aes(x = x,
                 y = y,
                 color = cu_type_color),
             size = 0.1) +
  scale_color_identity() +
  theme_classic(base_size = 5)

ggsave("cu_lowcat_tsne_cell_labels.pdf",
       cu_only_tsne,
       width = 2, height = 2,
       useDingbats = FALSE)
