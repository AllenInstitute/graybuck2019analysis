library(dplyr)
library(ggplot2)
options(stringsAsFactors = F)

# Common directory
cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/"

# Load metadata
samples <- read.csv(file.path(cdir,"all_samples.csv"), row.names = 1)

frag_cutoff <- 1e4
ENCODE_cutoff <- 0.25
chrom_cutoff <- 0.1

cutoff_binaries <- data.frame(pass_frag = as.numeric(samples$unique_fragments < frag_cutoff),
                              pass_ENCODE = as.numeric(samples$ENCODE_frac < ENCODE_cutoff),
                              pass_chrom = as.numeric(samples$frac_gt_250bp < chrom_cutoff))

cutoff_stats <- cutoff_binaries %>%
  group_by(pass_frag, pass_ENCODE, pass_chrom) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(-n) %>%
  mutate(xpos = 1:n()) %>%
  mutate(ypos_frag = ifelse(pass_frag, -250, NA),
         ypos_encode = ifelse(pass_ENCODE, -500, NA),
         ypos_chrom = ifelse(pass_chrom, -750, NA)) %>%
  rowwise() %>%
  mutate(ymin_line = min(ypos_frag, ypos_encode, ypos_chrom, na.rm = TRUE),
         ymax_line = max(ypos_frag, ypos_encode, ypos_chrom, na.rm = TRUE)) %>%
  mutate(ymin_line = ifelse(is.finite(ymin_line),
                            ymin_line, NA),
         ymax_line = ifelse(is.finite(ymax_line),
                            ymax_line, NA))

cutoff_sums <- data.frame(cutoff_label = c("<= 1e4 Unique Fragments",
                                           "<= 0.25 ENCODE DNase-Seq Overlap",
                                           "<= 0.1 fragments over 250bp"),
                          cutoff_total = c(sum(cutoff_binaries$pass_frag),
                                           sum(cutoff_binaries$pass_ENCODE),
                                           sum(cutoff_binaries$pass_chrom))) %>%
  mutate(cutoff_frac = cutoff_total/nrow(cutoff_binaries)) %>%
  mutate(ypos = c(-250, -500, -750))

upset_plot <- ggplot() +
  geom_rect(data = cutoff_stats,
            aes(xmin = xpos - 0.4,
                xmax = xpos + 0.4,
                ymin = 0,
                ymax = n)) +
  geom_rect(data = cutoff_sums,
            aes(xmin = -4 * cutoff_frac,
                xmax = 0,
                ymin = ypos - 100,
                ymax = ypos + 100)) +
  geom_text(data = cutoff_sums,
            aes(x = -2,
                y = ypos,
                label = cutoff_label),
            hjust = 1,
            vjust = 0.3) +
  geom_segment(data = cutoff_stats,
               aes(x = xpos,
                   xend = xpos,
                   y = ymin_line,
                   yend = ymax_line)) +
  geom_point(data = cutoff_stats,
             aes(x = xpos,
                 y = ypos_frag)) +
  geom_point(data = cutoff_stats,
             aes(x = xpos,
                 y = ypos_encode)) +
  geom_point(data = cutoff_stats,
             aes(x = xpos,
                 y = ypos_chrom)) +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_bw()

ggsave("fail_qc_upset.pdf",
       upset_plot,
       width = 4,
       height = 3,
       useDingbats = FALSE)
