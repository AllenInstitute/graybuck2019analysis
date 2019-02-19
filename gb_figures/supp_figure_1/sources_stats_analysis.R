library(dplyr)
library(ggplot2)
library(ggbeeswarm)
options(stringsAsFactors = F)

samples <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/all_samples.csv")


# Stats for driver lines
driver_samples <- samples %>%
  mutate(full_genotype = ifelse(is.na(inj_target),
                                full_genotype,
                                paste0("Z_Retrograde_",inj_target))) %>%
  mutate(full_genotype = ifelse(full_genotype == "Z_Retrograde_Retro_Orbital",
                                "ZZ_Retroorbital_mscRE4",
                                full_genotype)) %>%
  mutate(total_fragments = total_reads/2)

driver_xpos <- driver_samples %>%
  select(full_genotype) %>%
  unique() %>%
  arrange(full_genotype) %>%
  mutate(xpos = 1:n())

driver_samples <- driver_samples %>%
  left_join(driver_xpos)

# beeswarm for total_fragments
driver_total_stats <- driver_samples %>%
  group_by(full_genotype) %>%
  summarise(xpos = xpos[1],
            n = n(),
            median = median(total_fragments),
            q25 = quantile(total_fragments)[2],
            q75 = quantile(total_fragments)[4])

guide_rects <- driver_xpos %>%
  mutate(xmin = xpos - 0.5,
         xmax = xpos + 0.5,
         color = rep(c("#E6E7E8","#FFFFFF"), ceiling(n()/2))[1:n()],
         ymin = 0,
         ymax = max(driver_samples$total_fragments))

total_reads_jitter <- ggplot() +
  # Shaded guide rectangles
  geom_rect(data = guide_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = log10(ymin + 1),
                ymax = log10(ymax + 1),
                fill = color)) +
  # Jittered points for total_fragments
  geom_quasirandom(data = driver_samples,
                   aes(x = xpos,
                       y = log10(total_fragments + 1)),
                   color = "skyblue",
                   size = 0.1,
                   pch = 20) +
  # 1e4 cutoff line
  geom_hline(data = data.frame(y = 4),
             aes(yintercept = y),
             size = 0.5,
             linetype = "dashed") +
  # q25/q75 error bars
  geom_errorbar(data = driver_total_stats,
                aes(x = xpos,
                    ymin = log10(q25 + 1),
                    ymax = log10(q75 + 1)),
                size = 0.2,
                width = 0.4) +
  # Median points
  geom_point(data = driver_total_stats,
             aes(x = xpos,
                 y = log10(median + 1)),
             color = "red",
             size = 0.5) +
  # N cells
  geom_text(data = driver_total_stats,
            aes(x = xpos,
                y = 0.2,
                label = n),
            size = 1) +
  # Median values
  geom_text(data = driver_total_stats,
            aes(x = xpos,
                y = 0.75,
                label = round(log10(median + 1),1)),
            size = 1) +
  scale_fill_identity() +
  scale_x_reverse(breaks = driver_xpos$xpos,
                     labels = driver_xpos$full_genotype) +
  theme_classic(4) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.3)) +
  coord_flip()

ggsave("total_reads_jitter.pdf",
       total_reads_jitter,
       width = 2.5, height = 4,
       useDingbats = F)

# beeswarm for unique_fragments

driver_unique_stats <- driver_samples %>%
  group_by(full_genotype) %>%
  summarise(xpos = xpos[1],
            n = n(),
            median = median(unique_fragments),
            q25 = quantile(unique_fragments)[2],
            q75 = quantile(unique_fragments)[4])

guide_rects <- driver_xpos %>%
  mutate(xmin = xpos - 0.5,
         xmax = xpos + 0.5,
         color = rep(c("#E6E7E8","#FFFFFF"), ceiling(n()/2))[1:n()],
         ymin = 0,
         ymax = max(driver_samples$unique_fragments))


unique_fragments_jitter <- ggplot() +
  # Shaded guide rectangles
  geom_rect(data = guide_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = log10(ymin + 1),
                ymax = log10(ymax + 1),
                fill = color)) +
  # Jittered points for unique_fragments
  geom_quasirandom(data = driver_samples,
                   aes(x = xpos,
                       y = log10(unique_fragments + 1)),
                   color = "skyblue",
                   size = 0.1,
                   pch = 20) +
  # 1e4 cutoff line
  geom_hline(data = data.frame(y = 4),
             aes(yintercept = y),
             size = 0.5,
             linetype = "dashed") +
  # q25/q75 error bars
  geom_errorbar(data = driver_unique_stats,
                aes(x = xpos,
                    ymin = log10(q25 + 1),
                    ymax = log10(q75 + 1)),
                size = 0.2,
                width = 0.4) +
  # Median points
  geom_point(data = driver_unique_stats,
             aes(x = xpos,
                 y = log10(median + 1)),
             color = "red",
             size = 0.5) +
  # N cells
  geom_text(data = driver_unique_stats,
            aes(x = xpos,
                y = 0.2,
                label = n),
            size = 1) +
  # Median values
  geom_text(data = driver_unique_stats,
            aes(x = xpos,
                y = 0.75,
                label = round(log10(median + 1),1)),
            size = 1) +
  scale_fill_identity() +
  scale_x_reverse(breaks = driver_xpos$xpos,
                  labels = driver_xpos$full_genotype) +
  theme_classic(4) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.3)) +
  coord_flip()

ggsave("unique_fragments_jitter.pdf",
       unique_fragments_jitter,
       width = 2.5, height = 4,
       useDingbats = F)


# beeswarm for ENCODE_frac
driver_ENCODE_stats <- driver_samples %>%
  group_by(full_genotype) %>%
  summarise(xpos = xpos[1],
            n = n(),
            median = median(ENCODE_frac),
            q25 = quantile(ENCODE_frac)[2],
            q75 = quantile(ENCODE_frac)[4])

guide_rects <- driver_xpos %>%
  mutate(xmin = xpos - 0.5,
         xmax = xpos + 0.5,
         color = rep(c("#E6E7E8","#FFFFFF"), ceiling(n()/2))[1:n()],
         ymin = 0,
         ymax = max(driver_samples$ENCODE_frac))

ENCODE_frac_jitter <- ggplot() +
  # Shaded guide rectangles
  geom_rect(data = guide_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = color)) +
  # Jittered points for ENCODE_frac
  geom_quasirandom(data = driver_samples,
                   aes(x = xpos,
                       y = ENCODE_frac),
                   color = "skyblue",
                   size = 0.1,
                   pch = 20) +
  # 1e4 cutoff line
  geom_hline(data = data.frame(y = 0.25),
             aes(yintercept = y),
             size = 0.5,
             linetype = "dashed") +
  # q25/q75 error bars
  geom_errorbar(data = driver_ENCODE_stats,
                aes(x = xpos,
                    ymin = q25,
                    ymax = q75),
                size = 0.2,
                width = 0.4) +
  # Median points
  geom_point(data = driver_ENCODE_stats,
             aes(x = xpos,
                 y = median),
             color = "red",
             size = 0.5) +
  # N cells
  # geom_text(data = driver_ENCODE_stats,
  #           aes(x = xpos,
  #               y = 0.2,
  #               label = n),
  #           size = 1) +
  # Median values
  geom_text(data = driver_ENCODE_stats,
            aes(x = xpos,
                y = 0.05,
                label = round(median,2)),
            size = 1) +
  scale_fill_identity() +
  scale_x_reverse(breaks = driver_xpos$xpos,
                  labels = driver_xpos$full_genotype) +
  scale_y_continuous(limits = c(0, 0.75),
                     breaks = c(0, 0.25, 0.5, 0.75)) +
  theme_classic(4) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.3)) +
  coord_flip()

ggsave("ENCODE_fraction_jitter.pdf",
       ENCODE_frac_jitter,
       width = 2.5, height = 4,
       useDingbats = F)


# beeswarm for frac_gt_250bp
driver_chromatin_stats <- driver_samples %>%
  group_by(full_genotype) %>%
  summarise(xpos = xpos[1],
            n = n(),
            median = median(frac_gt_250bp),
            q25 = quantile(frac_gt_250bp)[2],
            q75 = quantile(frac_gt_250bp)[4])

guide_rects <- driver_xpos %>%
  mutate(xmin = xpos - 0.5,
         xmax = xpos + 0.5,
         color = rep(c("#E6E7E8","#FFFFFF"), ceiling(n()/2))[1:n()],
         ymin = 0,
         ymax = max(driver_samples$frac_gt_250bp))

frac_gt_250bp_jitter <- ggplot() +
  # Shaded guide rectangles
  geom_rect(data = guide_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = color)) +
  # Jittered points for frac_gt_250bp
  geom_quasirandom(data = driver_samples,
                   aes(x = xpos,
                       y = frac_gt_250bp),
                   color = "skyblue",
                   size = 0.1,
                   pch = 20) +
  # 1e4 cutoff line
  geom_hline(data = data.frame(y = 0.1),
             aes(yintercept = y),
             size = 0.5,
             linetype = "dashed") +
  # q25/q75 error bars
  geom_errorbar(data = driver_chromatin_stats,
                aes(x = xpos,
                    ymin = q25,
                    ymax = q75),
                size = 0.2,
                width = 0.4) +
  # Median points
  geom_point(data = driver_chromatin_stats,
             aes(x = xpos,
                 y = median),
             color = "red",
             size = 0.5) +
  # N cells
  # geom_text(data = driver_chromatin_stats,
  #           aes(x = xpos,
  #               y = 0.2,
  #               label = n),
  #           size = 1) +
  # Median values
  geom_text(data = driver_chromatin_stats,
            aes(x = xpos,
                y = 0.6,
                label = round(median,2)),
            size = 1) +
  scale_fill_identity() +
  scale_x_reverse(breaks = driver_xpos$xpos,
                  labels = driver_xpos$full_genotype) +
  scale_y_continuous(limits = c(0, 0.75),
                     breaks = c(0, 0.1, 0.25, 0.5, 0.75)) +
  theme_classic(4) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.3)) +
  coord_flip()

ggsave("chromatin_jitter.pdf",
       frac_gt_250bp_jitter,
       width = 2.5, height = 4,
       useDingbats = F)
