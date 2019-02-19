library(ggplot2)
library(dplyr)
library(scrattch.io)
options(stringsAsFactors)

# Analysis directory
adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-07-10_paper_analysis/"

# Common directory
cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/"

samples <- read.csv(file.path(cdir,"all_samples.csv"), row.names = 1)

# Load data frame with tSNE coordinates and phenograph clusters
mapping <- read.csv(file.path(adir,"2e4_tss_phenograph_scatac_visp_marker_correlation.csv"), row.names =1) %>%
  left_join(samples) %>%
  mutate(full_genotype = ifelse(is.na(inj_target),
                                full_genotype,
                                inj_target)) %>%
  annotate_cat(full_genotype)


foreground_points <- mapping %>%
  filter(full_genotype_label %in% c("Rbp4-Cre_KL100/wt;Ai14(RCL-tdT)/wt",
                             "LP",
                             "SC"))

background_points <- mapping %>%
  filter(!sample_id %in% foreground_points$sample_id)

# tSNE by pred_subclass
rbp4_tsne <- ggplot() +
  geom_point(data = background_points,
             aes(x = x,
                 y = y,
                 color = pred_subclass_color),
             size = 0.2,
             alpha = 0.5,
             pch = 16) +
  geom_point(data = foreground_points,
             aes(x = x,
                 y = y,
                 color = full_genotype_color),
             size = 0.4) +
  scale_color_identity() +
  scale_x_continuous("Lim1") +
  scale_y_continuous("Lim2")  +
  theme_bw(6) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.1,
                                        color = "#BEBEBE"),
        panel.border = element_blank(),
        axis.ticks = element_blank())

ggsave("rbp4_tsne.pdf",
       rbp4_tsne,
       width = 2.25,
       height = 2.25,
       useDingbats = F)

