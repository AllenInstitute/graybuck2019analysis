library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(lowcat)
library(cowplot)
options(stringsAsFactors = FALSE)

yshift_overlaps <- function(df, yinit = 1, yshift = -0.05, pad = 100) {
  df <- df %>%
    arrange(start, end) %>%
    mutate(y = yinit)
  
  for(i in 2:nrow(df)) {
    pre_df <- df[1:(i -1),]
    
    ol_df <- pre_df %>%
      filter(end + pad > df$start[i])
    if(nrow(ol_df) > 0) {
      # check for openings below
      ymax_ol <- max(ol_df$y)
      check_y <- seq(yinit, ymax_ol, by = yshift)
      open_y <- setdiff(check_y, ol_df$y)
      if(length(open_y) > 0) {
        df$y[i] <- open_y[1]
      } else {
        df$y[i] <- ymax_ol + yshift
      }
    }
    
  }
  df
}

build_process_plot <- function(fragment_GR,
                               target_region,
                               expand_size = 2e4,
                               frag_width = 1e3,
                               yshift = 0.05,
                               basecolor = "#808080",
                               samplecolor = "dodgerblue") {
  
  target_GR <- ucsc_loc_to_gr(target_region)
  
  exp_target_GR <- resize(target_GR, width = width(target_GR) + expand_size, fix = "center")
  
  fragment_GR_target <- subsetByOverlaps(fragment_GR, exp_target_GR)
  
  fragment_GR_target_exp <- resize(fragment_GR_target, width = frag_width, fix = "center")
  
  fragment_GR_target_exp_red <- reduce(fragment_GR_target_exp)
  
  fragment_GR_df <- as.data.frame(fragment_GR_target) %>% 
    yshift_overlaps(yinit = 3,
                    yshift = yshift)
  
  fragment_GR_df_exp <- as.data.frame(fragment_GR_target_exp) %>%
    yshift_overlaps(yinit = 2,
                    yshift = yshift)
  
  fragment_GR_df_exp_red <- as.data.frame(fragment_GR_target_exp_red) %>%
    mutate(y = 1)
  
  hlines <- data.frame(y_int = 1:3 - yshift)
  
  ggplot() +
    geom_hline(data = hlines,
               aes(yintercept = y_int),
               size = 0.5,
               color = basecolor) +
    geom_segment(data = fragment_GR_df,
                 aes(x = start, xend = end,
                     y = y, yend = y),
                 size = 1,
                 color = samplecolor) +
    geom_segment(data = fragment_GR_df_exp,
                 aes(x = start, xend = end,
                     y = y, yend = y),
                 size = 1,
                 color = samplecolor) +
    geom_segment(data = fragment_GR_df_exp_red,
                 aes(x = start, xend = end,
                     y = y, yend = y),
                 size = 1,
                 color = samplecolor) +
    scale_y_continuous("",limits = c(0.8, 3.5)) +
    scale_color_identity() +
    scale_x_continuous(seqnames(exp_target_GR),
                       limits = c(start(exp_target_GR), end(exp_target_GR))) +
    theme_bw(5) +
    theme(axis.ticks = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank())
}

build_tss_plot <- function(fragment_GR,
                           target_region,
                           target_strand = "+",
                           expand_size = 2e4,
                           yshift = 0.05,
                           basecolor = "#808080",
                           samplecolor = "dodgerblue") {
  
  target_GR <- ucsc_loc_to_gr(target_region)
  
  exp_target_GR <- resize(target_GR, width = width(target_GR) + expand_size, fix = "center")
  
  fragment_GR_target <- subsetByOverlaps(fragment_GR, exp_target_GR)
  
  fragment_GR_df <- as.data.frame(fragment_GR_target) %>% 
    yshift_overlaps(yinit = 1,
                    yshift = yshift)
  
  hlines <- data.frame(y_int = 1 - yshift)
  
  tss_pos <- GenomicRanges::start(target_GR)
  
  tss_seg <- data.frame(x = tss_pos,
                        xend = tss_pos,
                        y = 1 - yshift,
                        yend = 1.4)
  
  if(target_strand == "+") {
    tss_arrow <- data.frame(x = tss_pos,
                            xend = tss_pos + expand_size/5,
                            y = 1.4,
                            yend = 1.4)
  } else {
    tss_arrow <- data.frame(x = tss_pos - expand_size/5,
                            xend = tss_pos,
                            y = 1.4,
                            yend = 1.4)
  }

  
  ggplot() +
    geom_hline(data = hlines,
               aes(yintercept = y_int),
               size = 0.5,
               color = basecolor) +
    geom_segment(data = tss_seg,
                 aes(x = x, xend = xend,
                     y = y, yend = yend),
                 size = 0.55, 
                 color = basecolor,
                 lineend = "square") + 
    geom_segment(data = tss_arrow,
                 aes(x = x, xend = xend,
                     y = y, yend = yend),
                 size = 0.5,
                 color = basecolor,
                 arrow = arrow(type = "closed", length = unit(0.05,"inches"))) +
    geom_segment(data = fragment_GR_df,
                 aes(x = start, xend = end,
                     y = y, yend = y),
                 size = 0.5,
                 color = samplecolor) +
    scale_y_continuous("",limits = c(0.8, 1.5)) +
    scale_color_identity() +
    scale_x_continuous(seqnames(exp_target_GR),
                       limits = c(start(exp_target_GR), end(exp_target_GR))) +
    theme_bw(5) +
    theme(axis.ticks = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank())
}


# Process Plots

cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/"

samples <- read.csv(file.path(cdir, "all_samples.csv"))
qc_samples <- samples %>%
  filter(ENCODE_frac > 0.25,
         frac_gt_250bp > 0.1,
         unique_fragments > 1e4)

load(file.path(cdir,"bam_fragments.rda"))
bam_fragments <- bam_fragments[qc_samples$sample_id]

snap25_region <- "chr2:136,713,450-136,782,428"

p1 <- build_process_plot(bam_fragments[[1]],
                   snap25_region,
                   expand = 5e4,
                   frag_width = 5e3,
                   yshift = 0.15,
                   samplecolor = "dodgerblue")

p2 <- build_process_plot(bam_fragments[[2]],
                   snap25_region,
                   expand = 5e4,
                   frag_width = 5e3,
                   yshift = 0.15,
                   samplecolor = "orangered")

p3 <- build_process_plot(bam_fragments[[3]],
                   snap25_region,
                   expand = 5e4,
                   frag_width = 5e3,
                   yshift = 0.15,
                   samplecolor = "mediumorchid3")

all_plots <- plot_grid(p1, p2, p3,
                       ncol = 1)

ggsave("read_processing_plots.pdf",
       all_plots,
       width = 6, height = 4,
       useDingbats = F)

hspa8_region <- "chr9:40,801,273-40,805,199"

build_process_plot(bam_fragments[[3]],
                   hspa8_region)

# Cluster TSS plots

# Phenograph clustering
adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-07-10_paper_analysis/"

pmc <- read.csv(file.path(adir, "2e4_tss_phenograph_scatac_visp_correlation.csv"))
pmc_clusters <- pmc %>%
  group_by(pg_clust,pg_clust_color, pred_cluster_id, pred_cluster_label) %>%
  summarise(n_samples = n()) %>%
  arrange(pred_cluster_id)

snap25_tss <- "chr2:136,713,450-136,713,460"
slc17a7_tss <- "chr7:45,163,921-45,163,922"

# 86: Lamp5 Lhx6
# 19: L5 PT VISp C1ql2 Cdh13
# 58: Oligo Serpinb1a_2
# 61: L6b

pmc_keep <- pmc %>%
  filter(pg_clust %in% c(86, 61, 58))

pmc_fragments <- combine_group_GRanges(bam_fragments,
                                       pmc_keep,
                                       "pg_clust")

lamp5_snap25 <- build_tss_plot(fragment_GR = pmc_fragments[["86"]],
               target_region = snap25_tss,
               expand_size = 2.2e4,
               yshift = 0.05,
               basecolor = "#808080",
               samplecolor = pmc_clusters$pg_clust_color[pmc_clusters$pg_clust == "86"])

l6b_snap25 <- build_tss_plot(fragment_GR = pmc_fragments[["61"]],
                              target_region = snap25_tss,
                              expand_size = 2.2e4,
                              yshift = 0.05,
                              basecolor = "#808080",
                              samplecolor = pmc_clusters$pg_clust_color[pmc_clusters$pg_clust == "61"])

oligo_snap25 <- build_tss_plot(fragment_GR = pmc_fragments[["58"]],
               target_region = snap25_tss,
               expand_size = 2.2e4,
               yshift = 0.05,
               basecolor = "#808080",
               samplecolor = pmc_clusters$pg_clust_color[pmc_clusters$pg_clust == "58"])

lamp5_slc17a7 <- build_tss_plot(fragment_GR = pmc_fragments[["86"]],
                               target_region = slc17a7_tss,
                               expand_size = 2.2e4,
                               yshift = 0.05,
                               basecolor = "#808080",
                               samplecolor = pmc_clusters$pg_clust_color[pmc_clusters$pg_clust == "86"])

l6b_slc17a7 <- build_tss_plot(fragment_GR = pmc_fragments[["61"]],
                             target_region = slc17a7_tss,
                             expand_size = 2.2e4,
                             yshift = 0.05,
                             basecolor = "#808080",
                             samplecolor = pmc_clusters$pg_clust_color[pmc_clusters$pg_clust == "61"])

oligo_slc17a7 <- build_tss_plot(fragment_GR = pmc_fragments[["58"]],
                               target_region = slc17a7_tss,
                               expand_size = 2.2e4,
                               yshift = 0.05,
                               basecolor = "#808080",
                               samplecolor = pmc_clusters$pg_clust_color[pmc_clusters$pg_clust == "58"])

all_tss_plots <- plot_grid(lamp5_snap25, l6b_snap25, oligo_snap25,
                           lamp5_slc17a7, l6b_slc17a7, oligo_slc17a7,
                           ncol = 3)

ggsave("tss_cluster_counts_processing_plots.pdf",
       all_tss_plots,
       width = 6, height = 2,
       useDingbats = F)
