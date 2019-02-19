library(ggplot2)
library(scrattch.vis)
options(stringsAsFactors = F)

load("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-06-26_aws_parameter_testing/jaccard/f1e4_e25_c10_ds1de4_x1e3_jaccard.rda")

cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/"

jaccard_subset <- jaccard_results %>%
  filter(s1 < 201 & s2 < 201)

jaccard_subset_matrix <- jaccard_matrix[1:200,1:200]

# save the jaccard matrix subset to CSV
write.csv(jaccard_subset_matrix,
          "Fig_2_Jaccard_distance_matrix_subset.csv")

jaccard_subset_hc <- hclust(as.dist(jaccard_subset_matrix), method = "ward.D")

plot_data <- jaccard_subset %>%
  mutate(color = values_to_colors(jaccard_distance,
                                  minval = 0.92,
                                  maxval = 0.99,
                                  colorset = c("black","white"))) %>%
  rowwise() %>%
  mutate(s1_pos = which(jaccard_subset_hc$order == s1),
         s2_pos = which(jaccard_subset_hc$order == s2))

distance_plot <- ggplot() +
  geom_tile(data = plot_data,
            aes(x = s1_pos,
                y = -s2_pos,
                fill = color)) +
  geom_tile(data = plot_data,
            aes(x = s2_pos,
                y = -s1_pos,
                fill = color)) +
  scale_fill_identity() +
  theme_void()

ggsave("distance_heatmap.png",
       distance_plot,
       width = 1.5, height = 1.5)
