library(ggplot2)
library(dplyr)
library(scrattch.vis)
library(scrattch.io)
options(stringsAsFactors = FALSE)

# Analysis directory
adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-07-10_paper_analysis/"
cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/"

samples <- read.csv(file.path(cdir,"all_samples.csv"), row.names = 1)
source_cats <- read.csv("source_category_colors.csv", comment.char = "|")

samples <- left_join(samples, source_cats)

# Load data frame with tSNE coordinates
load(file.path(adir,"f1e4_e25_c10_ds1e4_x1e3_jaccard.rda"))

plot_data <- tsne_df %>%
  left_join(samples) %>%
  #left_join(source_cats) %>%
  annotate_cat("full_genotype")

fg_colors <- plot_data %>%
  select(full_genotype_label, full_genotype_color) %>%
  unique()


# tSNE by source category
source_tsne <- ggplot() +
  geom_point(data = plot_data,
             aes(x = x,
                 y = y,
                 color = source_category_color),
             size = 0.05) +
  scale_color_identity() +
  scale_x_continuous("Lim1") +
  scale_y_continuous("Lim2")  +
  theme_bw(6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.1,
                                        color = "#BEBEBE"),
        panel.border = element_blank(),
        axis.ticks = element_blank())

ggsave("source_tsne.pdf",
       source_tsne,
       width = 2.25,
       height = 2.25,
       useDingbats = F)

# Centroid classifier results
mapping <- read.csv(file.path(adir,"2e4_tss_phenograph_scatac_visp_marker_correlation.csv"), row.names = 1) %>%
  select(-x, -y)

plot_data <- left_join(plot_data, mapping)

# tSNE by subclass
subclass_tsne <- ggplot() +
  geom_point(data = plot_data,
             aes(x = x,
                 y = y,
                 color = pred_subclass_color),
             size = 0.05) +
  scale_color_identity() +
  scale_x_continuous("Lim1") +
  scale_y_continuous("Lim2")  +
  theme_bw(6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.1,
                                        color = "#BEBEBE"),
        panel.border = element_blank(),
        axis.ticks = element_blank())

ggsave("subclass_tsne.pdf",
       subclass_tsne,
       width = 2.25,
       height = 2.25,
       useDingbats = F)

# tSNE by cluster
cluster_tsne <- ggplot() +
  geom_point(data = plot_data,
             aes(x = x,
                 y = y,
                 color = pred_cluster_color),
             size = 0.05) +
  scale_color_identity() +
  scale_x_continuous("Lim1") +
  scale_y_continuous("Lim2")  +
  theme_bw(6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.1,
                                        color = "#BEBEBE"),
        panel.border = element_blank(),
        axis.ticks = element_blank())

ggsave("cluster_tsne.pdf",
       cluster_tsne,
       width = 2.25,
       height = 2.25,
       useDingbats = F)

# tSNE by cluster with labels
pred_cluster_plot_labels <- plot_data %>%
  select(x, y, pred_cluster_label) %>%
  group_by(pred_cluster_label) %>%
  summarise(x = mean(x),
            y = mean(y)) %>%
  ungroup()

cluster_tsne_labels <- ggplot() +
  geom_point(data = plot_data,
             aes(x = x,
                 y = y,
                 color = pred_cluster_color),
             size = 0.05) +
  geom_text(data = pred_cluster_plot_labels,
            aes(x = x,
                y = y,
                label = pred_cluster_label)) +
  scale_color_identity() +
  scale_x_continuous("Lim1") +
  scale_y_continuous("Lim2")  +
  theme_bw(6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.1,
                                        color = "#BEBEBE"),
        panel.border = element_blank(),
        axis.ticks = element_blank())

ggsave("cluster_tsne_labels.pdf",
       cluster_tsne_labels,
       width = 2.25,
       height = 2.25,
       useDingbats = F)

# tSNE by phenograph cluster

pg_clust_tsne <- ggplot() +
  geom_point(data = plot_data,
             aes(x = x,
                 y = y,
                 color = pg_clust_color),
             size = 0.05) +
  scale_color_identity() +
  scale_x_continuous("Lim1") +
  scale_y_continuous("Lim2")  +
  theme_bw(6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.1,
                                        color = "#BEBEBE"),
        panel.border = element_blank(),
        axis.ticks = element_blank())

ggsave("pg_clust_tsne.pdf",
       pg_clust_tsne,
       width = 2.25,
       height = 2.25,
       useDingbats = F)


# Save data for tSNE plots

req_plot_data <- plot_data %>%
  select(sample_id, x, y, 
         source_category, source_category_color,
         pg_clust, pg_clust_color,
         pred_subclass_label, pred_subclass_cor, pred_subclass_color,
         pred_cluster_label, pred_cluster_cor, pred_cluster_color)

write.csv(req_plot_data,
          "Fig_2_tSNE_plot_data.csv",
          row.names = FALSE)
