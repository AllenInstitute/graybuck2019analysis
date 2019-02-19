library(ggplot2)
library(dplyr)
library(reshape2)
library(scrattch.io)
library(scrattch.vis)
options(stringsAsFactors = F)

adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-06-19_1kb_fragments/"

load(file.path(adir, "pg_cluster_cor_results.rda"))

load(file.path(adir, "pg_tsne_plot_annotations.rda"))

pg_cluster_anno <- pg_tsne_plot_data %>%
  select(pg_cluster_id, pg_cluster_color, cluster_id, cluster_cor) %>%
  unique() %>%
  arrange(cluster_id) %>%
  mutate(ypos = n():1) %>%
  select(-cluster_id)

anno <- read_tome_anno("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/tomes/facs/mouse_V1_ALM_20180520/transcrip.tome")
anno <- anno %>%
  # Remove outlier clusters
  filter(cluster_id %in% 1:133) %>%
  # Filter ALM-only clusters
  group_by(cluster_id) %>%
  mutate(alm_frac = sum(region_label == "ALM")/n()) %>%
  ungroup() %>%
  filter(alm_frac < 0.9) %>%
  # Only keep VISp cells
  filter(region_label == "VISp")

cluster_anno <- anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique() %>%
  arrange(cluster_id) %>%
  mutate(xpos = 1:n())

max_pg_clusters <-  pg_tsne_plot_data %>%
  select(pg_cluster_id, pg_cluster_color, cluster_id, cluster_cor) %>%
  unique() %>%
  arrange(cluster_id) %>%
  mutate(ypos = n():1) %>%
  left_join(cluster_anno)

scaled_all_pg_cluster_cor_df <- t(apply(all_pg_cluster_cor_df,
                                      1,
                                      function(x) {
                                        x/max(x)
                                      }))

all_cor_plot_data <- melt(scaled_all_pg_cluster_cor_df)
names(all_cor_plot_data) <- c("pg_cluster_id","cluster_label","cor")
all_cor_plot_data <- all_cor_plot_data %>%
  left_join(pg_cluster_anno) %>%
  left_join(cluster_anno)

colorset <- c("dodgerblue","gray80","orangered")

all_cor_plot_data <- all_cor_plot_data %>%
  mutate(fill = values_to_colors(cor + -1*min(cor), colorset = colorset))

ggplot() +
  geom_tile(data = all_cor_plot_data,
            aes(x = xpos,
                y = ypos,
                fill = fill))+
    geom_tile(data = max_pg_clusters,
              aes(x = xpos, y = ypos),
              fill = NA,
              color = "white",
              size = 0.1) +
  scale_fill_identity() +
  scale_x_discrete("scRNA-seq Cluster", expand = c(0,0)) +
  scale_y_discrete("scATAC-seq Cluster", expand = c(0,0))
