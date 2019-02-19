library(lowcat)
library(dplyr)
library(scrattch.io)
library(scrattch.vis)
library(GenomicRanges)
library(ggplot2)
options(stringsAsFactors = F)

# Analysis directory
adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-07-10_paper_analysis/"
# Common directory
cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common"

# Load bam fragments
load(file.path(cdir,"bam_fragments.rda"))
mapping <- read.csv(file.path(adir,"2e4_tss_phenograph_scatac_visp_marker_correlation.csv"))

# Set up custom subclass + distinct types for plotting
group_anno <- read.csv("cluster_grouping_for_tracks.csv", header = TRUE) %>%
  select(-pred_cluster_label, -pred_cluster_color)

mapping <- mapping %>%
  left_join(group_anno) %>%
  arrange(-group_id)

bam_fragments <- bam_fragments[mapping$sample_id]

mscre13_loc <- "chr12:9,487,900-9,488,400"
Osr1_loc <- "chr12:9,574,442-9,581,500"

group_anno <- group_anno %>%
  select(group_id, group_label, group_color) %>%
  unique() %>%
  filter(group_id %in% mapping$group_id) %>%
  arrange(-group_id)

group_colors <- group_anno$group_color
names(group_colors) <- group_anno$group_label


broad_plot <- build_pile_plot(bam_fragments,
                              Osr1_loc,
                              highlight_loc = mscre13_loc,
                              padding = c(1e5,2e4),
                              gr_groups = mapping$group_label,
                              group_colors = group_colors,
                              highlight_color = "#FFB517",
                              max_val = 5,
                              window_size = 100,
                              window_mode = "max")

ggsave("mscre13_broad_plot.pdf",
       broad_plot,
       width = 5.5,
       height = 3,
       useDingbats = F)

ggsave("mscre13_broad_plot.png",
       broad_plot,
       width = 5.5,
       height = 3)

