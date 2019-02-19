library(dplyr)
library(ggplot2)
library(ggExtra)

adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-06-07_all_mouse/"

samples <- read.csv(file.path(adir,"all_samples.csv"), row.names = 1)

samples <- samples %>%
  filter(ENCODE_frac < 1) %>%
  mutate(color = ifelse(ENCODE_frac > 0.25 & unique_fragments > 1e4,
                        "#009444",
                        "#000000"))

summaries <- samples %>%
  mutate(pass_encode = ENCODE_frac > 0.25) %>%
  mutate(pass_fragments = unique_fragments > 1e4) %>%
  group_by(pass_encode, pass_fragments) %>%
  summarise(n_cells = n(),
         perc_cells = round(n()/nrow(samples)*100,2)) %>%
  mutate(x = ifelse(pass_fragments, 6, 1.5),
         y = ifelse(pass_encode, 0.7, 0.1))

p <- ggplot()  +
  geom_point(data = samples,
             aes(x = log10(unique_fragments),
                 y = ENCODE_frac,
                 color = color),
             size = 0.5) +
  geom_segment(aes(x = 4, xend = 4,
                   y = 0, yend = 0.75),
            linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 6.5,
                   y = 0.25, yend = 0.25),
               linetype = "dashed") +
  geom_text(data = summaries,
            aes(x = x,
                y = y,
                label = paste0(perc_cells,"%")),
            hjust = 1) +
  geom_text(data = summaries,
            aes(x = x,
                y = y - 0.05,
                label = paste(n_cells,"cells")),
            hjust = 1) +
  theme_bw() +
  scale_color_identity() +
  scale_x_continuous("Log10(Unique Fragments)",
                     expand = c(0,0),
                     limits = c(-0.3, 6.5)) +
  scale_y_continuous("Fraction overlapping ENCODE peaks",
                     expand = c(0,0),
                     limits = c(-0.02, 0.75))

pm <- ggMarginal(p, 
                 type = "histogram",
                 fill = "#BEBEBE",
                 xparams = list(binwidth = 0.2),
                 yparams = list(binwidth = 0.025))

pm

ggsave("unique_vs_ENCODE.pdf",
       pm,
       width = 5,
       height = 5,
       useDingbats = F)
