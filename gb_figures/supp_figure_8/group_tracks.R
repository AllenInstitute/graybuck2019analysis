library(lowcat)
library(dplyr)
library(scrattch.io)
library(scrattch.vis)
library(GenomicRanges)
library(ggplot2)
library(purrr)
library(rtracklayer)
options(stringsAsFactors = F)
source("score_track_functions.R")

# Analysis directory
adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-07-10_paper_analysis/"
# Common directory
cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common"

# Load bam fragments
load(file.path(cdir,"bam_fragments.rda"))
mapping <- read.csv(file.path(adir,"2e4_tss_phenograph_scatac_visp_marker_correlation.csv"))

# Load peak locations
mscre_locs <- read.table("mscre_locs.txt", sep = "\t", header = TRUE)
mscre_bed <- map_dfr(mscre_locs$ucsc_loc, ucsc_loc_to_bed)
mscre_bed$name <- mscre_locs$id
mscre_bed$score <- 0
mscre_gr <- bed_to_GRanges(mscre_bed)

import_gr <- resize(mscre_gr, 5e4, fix = "center")

# Load PhastCons
ucsc_pc <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/mm10.60way.phastCons60wayEuarchontoGlire.bw"

pc_bw_import <- import(ucsc_pc, format = "bw", which = import_gr)

# Set up custom subclass + distinct types for plotting
group_anno <- read.csv("cluster_grouping_for_tracks.csv", header = TRUE) %>%
  select(-pred_cluster_label, -pred_cluster_color)

mapping <- mapping %>%
  left_join(group_anno) %>%
  arrange(-group_id) %>%
  select(sample_id, group_id, group_label, group_color)

bam_fragments <- bam_fragments[mapping$sample_id]

group_colors <- unique(mapping$group_color)
names(group_colors) <- unique(mapping$group_label)

# Generate plots

all_plots <- list()

for(x in 1:nrow(mscre_locs)) {
  broad_plot <- build_pile_plot(gr_list = bam_fragments,
                                ucsc_loc = mscre_locs$ucsc_loc[x],
                                highlight_loc = mscre_locs$ucsc_loc[x],
                                padding = c(2e4,2e4),
                                gr_groups = mapping$group_label,
                                group_colors = group_colors,
                                highlight_color = "#FFB517",
                                norm = "PM",
                                max_val = 4,
                                window_size = 20,
                                window_mode = "max")
  
  broad_plot <- add_score_track(broad_plot,
                                score_gr = pc_bw_import,
                                ucsc_loc = mscre_locs$ucsc_loc[x],
                                highlight_loc = mscre_locs$ucsc_loc[x],
                                track_ypos = length(group_colors) + 1,
                                padding = c(2e4,2e4),
                                highlight_color = "#FFB517",
                                norm = "max",
                                max_val = 1,
                                window_size = 20,
                                window_mode = "max")
  
  all_plots <- c(all_plots, list(broad_plot))
  
  ggsave(paste0(mscre_locs$id[x],"_tracks.pdf"),
         broad_plot,
         width = 5.5,
         height = 3,
         useDingbats = F)
  
  narrow_plot <- build_pile_plot(gr_list = bam_fragments,
                                 ucsc_loc = mscre_locs$ucsc_loc[x],
                                 highlight_loc = mscre_locs$ucsc_loc[x],
                                 padding = c(1e3,1e3),
                                 gr_groups = mapping$group_label,
                                 group_colors = group_colors,
                                 highlight_color = "#FFB517",
                                 norm = "none",
                                 max_val = NULL,
                                 window_size = 1,
                                 window_mode = "max")
  
  narrow_plot <- add_score_track(narrow_plot,
                                 score_gr = pc_bw_import,
                                 ucsc_loc = mscre_locs$ucsc_loc[x],
                                 highlight_loc = mscre_locs$ucsc_loc[x],
                                 track_ypos = length(group_colors) + 1,
                                 padding = c(1e3,1e3),
                                 highlight_color = "#FFB517",
                                 norm = "max",
                                 max_val = 1,
                                 window_size = 1,
                                 window_mode = "max")
  
  all_plots <- c(all_plots, list(narrow_plot))
  
  ggsave(paste0(mscre_locs$id[x],"_narrow_tracks.pdf"),
         narrow_plot,
         width = 2,
         height = 3,
         useDingbats = F)
  
}

library(cowplot)

joint_plot <- plot_grid(plotlist = all_plots,
                       ncol = 4,
                       rel_widths = c(2,1,2,1))

ggsave("all_track_plots.pdf",
       joint_plot,
       width = 7.5,
       height = 10)
