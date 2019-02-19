library(dplyr)
library(purrr)
library(ggplot2)
library(GenomicRanges)
library(matrixStats)
options(stringsAsFactors = F)

cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/"

load(file.path(cdir,"bam_fragments.rda"))
all_samples <- read.csv(file.path(cdir,"all_samples.csv"))

fragment_size_df <- function(frag_GR) {
  size_df <- as.data.frame(table(width(frag_GR)))
  names(size_df) <- c("width","n_fragments")
  size_df$width <- as.numeric(size_df$width)
  
  size_df <- size_df %>%
    mutate(frac_fragments = n_fragments/sum(n_fragments))
  size_df
}

fragment_size_matrix <- function(gr_list,
                                 min = 0,
                                 max = 1000) {
  size_mat <- matrix(0, nrow = length(min:max), ncol = length(gr_list))
  rownames(size_mat) <- min:max
  colnames(size_mat) <- names(gr_list)
  size_dfs <- map(gr_list, fragment_size_df)
  size_dfs <- map(size_dfs,
                  function(x) {
                    x %>%
                      filter(width >= min) %>%
                      filter(width <= max)
                  })
  for(i in 1:length(size_dfs)) {
    size_mat[size_dfs[[i]]$width, i] <- size_dfs[[i]]$frac_fragments
  }
  size_mat
}

pass_samples <- all_samples %>%
  filter(frac_gt_250bp > 0.1)

pass_size_mat <- fragment_size_matrix(bam_fragments[pass_samples$sample_id])

pass_df <- data.frame(size = as.numeric(rownames(pass_size_mat)),
                         med = rowMedians(pass_size_mat),
                         q25 = apply(pass_size_mat, 1, function(x) {
                           quantile(x)[2]
                         }),
                         q75 = apply(pass_size_mat, 1, function(x) {
                           quantile(x)[4]
                         }))

pass_plot <- ggplot() +
  geom_ribbon(data = pass_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "darkgreen",
              alpha = 0.4) +
  geom_line(data = pass_df,
            aes(x = size,
                y = med),
            color = "darkgreen",
            size = 0.1) +
  theme_bw(4)

ggsave("insert_size_pass.pdf",
       pass_plot,
       width = 1.5,
       height = 1)

fail_samples <- all_samples %>%
  filter(frac_gt_250bp <= 0.1)

fail_size_mat <- fragment_size_matrix(bam_fragments[fail_samples$sample_id])

fail_df <- data.frame(size = as.numeric(rownames(fail_size_mat)),
                      med = rowMedians(fail_size_mat),
                      q25 = apply(fail_size_mat, 1, function(x) {
                        quantile(x)[2]
                      }),
                      q75 = apply(fail_size_mat, 1, function(x) {
                        quantile(x)[4]
                      }))

fail_plot <- ggplot() +
  geom_ribbon(data = fail_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "gray80",
              alpha = 0.7) +
  geom_line(data = fail_df,
            aes(x = size,
                y = med),
            color = "black",
            size = 0.1) +
  theme_bw(4)

ggsave("insert_size_fail.pdf",
       fail_plot,
       width = 1.5,
       height = 1)
