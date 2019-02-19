library(scrattch.io)
library(scrattch.vis)
library(dplyr)
library(ggplot2)
library(lowcat)
library(reshape2)
library(viridisLite)
library(matrixStats)
options(stringsAsFactors = F)

# Load cluster medians
load("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/median_matrixes.rda")
# Load cluster annotations
# Load annotations
anno <- read_tome_anno("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/tomes/facs/mouse_V1_ALM_20180520/transcrip.tome")
# Remove outlier clusters
anno <- anno %>%
  filter(cluster_id %in% 1:133) %>%
  # Remove ALM-only clusters
  group_by(cluster_id) %>%
  mutate(alm_frac = sum(region_label == "ALM")/n()) %>%
  ungroup() %>%
  filter(alm_frac < 0.9) %>%
  filter(region_label == "VISp")

cluster_anno <- anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique() %>%
  arrange(cluster_id) %>%
  mutate(xpos = 1:n())

nn_labels <- anno %>%
  filter(class_label %in% c("Endothelial","Non-Neuronal")) %>%
  select(cluster_label) %>%
  unique() %>%
  unlist()

# Load accessibility scores per TSS
load("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-07-10_paper_analysis/2e4_tss_phenograph_scatac_visp_marker_correlation_matrices.rda")
# Load mapping results for colors
adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-07-10_paper_analysis/"

cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/"

pmc <- read.csv(file.path(adir, "2e4_tss_phenograph_scatac_visp_marker_correlation.csv"))
pmc_clusters <- pmc %>%
  group_by(pg_clust,pg_clust_color, pred_cluster_id, pred_cluster_label) %>%
  summarise(n_samples = n()) %>%
  arrange(pred_cluster_id)


cluster_med_mat <- cluster_med_mat[rownames(cluster_tss_counts),]

marker_df <- read.csv("markers_2018-05-25.csv")
selected_markers <- c(marker_df$Global, marker_df$Inh, marker_df$Ex)
selected_markers <- selected_markers[selected_markers != ""]
# Remove markers that aren't in both datasets
selected_markers <- selected_markers[selected_markers %in% rownames(cluster_med_mat) & selected_markers %in% rownames(subclass_tss_counts)]
# Remove markers that have 0 for medians
selected_markers <- selected_markers[rowMaxs(cluster_med_mat[selected_markers,]) > 0]
# Pick every 3rd gene - otherwise it looks noisy because the plot is so small.
selected_markers <- selected_markers[1:length(selected_markers) %% 3 == 0]

marker_med_mat <- cluster_med_mat[selected_markers,]

# Save the raw median matrix to csv
write.csv(marker_med_mat,
          "Fig_2_scRNAseq_marker_cluster_medians.csv")

gene_ypos <- data.frame(ypos = length(selected_markers):1,
                        gene = selected_markers)

# Normalize per-gene to make heatmap clearer
selected_med_mat <- marker_med_mat[, cluster_anno$cluster_label]
selected_med_mat <- selected_med_mat / rowMaxs(selected_med_mat)

cluster_med_melt <- melt(selected_med_mat)

cluster_med_plot_data <- cluster_med_melt %>%
  rename_("gene" = "Var1",
          "cluster_label" = "Var2") %>%
  left_join(cluster_anno) %>%
  left_join(gene_ypos) %>%
  mutate(fill = values_to_colors(log10(value + 1)))

exp_heatmap <- ggplot() +
  geom_tile(data = cluster_med_plot_data,
            aes(x = xpos,
                y = ypos,
                fill = fill)) +
  geom_rect(data = cluster_anno,
            aes(xmin = xpos - 0.5,
                xmax = xpos + 0.5,
                ymin = -4.5,
                ymax = 0.5,
                fill = cluster_color)) +
  scale_fill_identity() +
  theme_void()

exp_heatmap

ggsave("manual100_marker_expression_heatmap.pdf",
       exp_heatmap,
       width = 1.5, height = 1.5,
       useDingbats = F)

## Phenograph cluster tss counts

marker_tss_counts <- cluster_tss_counts[selected_markers,]

# Save the raw marker median matrix to CSV
write.csv(marker_tss_counts,
          "Fig_2_scATACseq_marker_cluster_counts.csv")

# Normalize per cluster to account for differences in depth between clusters
marker_tss_counts <- apply(marker_tss_counts, 2, scale)
rownames(marker_tss_counts) <- selected_markers

pg_xpos <- pmc_clusters %>%
  ungroup() %>%
  arrange(pred_cluster_id) %>%
  mutate(xpos = 1:n()) %>%
  select(pg_clust, xpos)

cluster_tss_melt <- melt(marker_tss_counts)

cluster_tss_counts_plot_data <- cluster_tss_melt %>%
  rename_("gene" = "Var1",
          "pg_clust" = "Var2") %>%
  left_join(pg_xpos) %>%
  left_join(gene_ypos) %>%
  mutate(fill = values_to_colors(value, colorset = viridis(10)))

pmc_clusters <- pmc_clusters %>%
  left_join(pg_xpos)
  

atac_heatmap <- ggplot() +
  geom_tile(data = cluster_tss_counts_plot_data,
            aes(x = xpos,
                y = ypos,
                fill = fill)) +
  geom_rect(data = pmc_clusters,
            aes(xmin = xpos - 0.5,
                xmax = xpos + 0.5,
                ymin = -4.5,
                ymax = 0.5,
                fill = pg_clust_color)) +
  scale_fill_identity() +
  theme_void()

atac_heatmap

ggsave("manual100_marker_atac_heatmap.pdf",
       atac_heatmap,
       width = 1.5, height = 1.5,
       useDingbats = F)
