library(dplyr)
library(ggplot2)
library(scrattch.io)
library(scrattch.vis)
library(purrr)
library(cowplot)
options(stringsAsFactors = F)

# Analysis directory
adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-07-10_paper_analysis/"
# Common directory
cdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/"

# Load metadata
samples <- read.csv(file.path(cdir,"all_samples.csv"), row.names = 1)

grouping <- read.csv("cluster_grouping_for_tracks.csv")

# Fix a couple of genotypes to match transcriptomic annotations
samples <- samples %>%
  mutate(full_genotype = ifelse(full_genotype == "Nkx2-CreERT2;Ai14(TAM at e17)",
                                "Nkx2-1-CreERT2/wt;Ai14(RCL-tdT)/wt",
                                full_genotype),
         full_genotype = ifelse(full_genotype == "Ntsr1-Cre/wt;Ai14(RCL-tdT)/wt",
                                 "Ntsr1-Cre_GN220/wt;Ai14(RCL-tdT)/wt",
                                 full_genotype),
         full_genotype = ifelse(full_genotype == "Slc17a8-IRES2-Cre;Ai14",
                                "Slc17a8-iCre/wt;Ai14(RCL-tdT)/wt",
                                full_genotype),
         full_genotype = ifelse(full_genotype == "Tac1-IRES-Cre-D/wt;Ai14(RCL-tdT)/wt",
                                "Tac1-IRES2-Cre/wt;Ai14(RCL-tdT)/wt",
                                full_genotype),
         full_genotype = ifelse(full_genotype == "Vipr2-IRES2-Cre;Slc32a1-T2A-FlpO;Ai65",
                                "Vipr2-IRES2-Cre/wt;Slc32a1-T2A-FlpO/wt;Ai65(RCFL-tdT)/wt",
                                full_genotype))

# Load mapping
mapping <- read.csv(file.path(adir,"2e4_tss_phenograph_scatac_visp_marker_correlation.csv"), row.names = 1)

mapping <- mapping %>%
  left_join(samples) %>%
  filter(is.na(inj_target)) %>%
  annotate_cat("full_genotype") %>%
  left_join(grouping)

genotype_pos <- mapping %>%
  select(full_genotype_label, full_genotype_id) %>%
  unique()
names(genotype_pos) <- c("genotype_label","genotype_pos")

genotypes <- genotype_pos$full_genotype_label

mapping_rect_data <- mapping %>%
  group_by(full_genotype_label) %>%
  mutate(genotype_n = n()) %>%
  group_by(full_genotype_id, full_genotype_label, full_genotype_color, genotype_n, group_id, group_label, group_color) %>%
  summarise(group_n = n()) %>%
  ungroup() %>%
  mutate(group_frac = group_n / genotype_n) %>%
  arrange(group_id) %>%
  group_by(full_genotype_label) %>%
  mutate(group_max = cumsum(group_frac)) %>%
  mutate(group_min = lag(group_max, default = 0)) %>%
  ungroup()

mapping_plot_data <- mapping_rect_data %>%
  select(full_genotype_label, group_min, group_max, group_color) %>%
  mutate(modality = "scATAC")

names(mapping_plot_data) <- c("genotype_label","group_min","group_max","group_color", "modality")

# Load transcriptomic annotations
tome <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/DataRelease_June2018/shinyapps.io_tome_files/Mouse_ALM_VISp/transcrip.tome"

anno <- read_tome_anno(tome)

anno <- anno %>%
  filter(brain_region_label == "VISp") %>%
  mutate(genotype_label = ifelse(genotype_label == "Rbp4-Cre_KL1/wt;Ai14(RCL-tdT)/wt",
                                 "Rbp4-Cre_KL100/wt;Ai14(RCL-tdT)/wt",
                                 genotype_label),
         genotype_label = ifelse(genotype_label == "Ntsr1-Cre_GN22/wt;Ai14(RCL-tdT)/wt",
                                 "Ntsr1-Cre_GN220/wt;Ai14(RCL-tdT)/wt",
                                 genotype_label))

matching_genotypes <- intersect(genotypes, anno$genotype_label)

# Removing Slc17a7, because it was collected as tdT(-)
matching_genotypes <- matching_genotypes[!grepl("Slc17a7",matching_genotypes)]

missing_genotypes <- setdiff(genotypes, anno$genotype_label)

txn_rect_data <- anno %>%
  filter(genotype_label %in% matching_genotypes) %>%
  group_by(genotype_label) %>%
  mutate(genotype_n = n()) %>%
  group_by(genotype_id, genotype_label, genotype_color, genotype_n, group_id, group_label, group_color) %>%
  summarise(group_n = n()) %>%
  ungroup() %>%
  mutate(group_frac = group_n / genotype_n) %>%
  arrange(group_id) %>%
  group_by(genotype_label) %>%
  mutate(group_max = cumsum(group_frac)) %>%
  mutate(group_min = lag(group_max, default = 0)) %>%
  ungroup()

txn_plot_data <- txn_rect_data %>%
  select(genotype_label, group_min, group_max, group_color) %>%
  mutate(modality = "scRNA")

modality_pos <- data.frame(ypos = 1:2,
                           modality = c("scATAC","scRNA"))

all_plot_data <- rbind(mapping_plot_data, txn_plot_data) %>%
  left_join(genotype_pos) %>%
  left_join(modality_pos)

group_bars <- ggplot(data = all_plot_data) +
  geom_rect(aes(xmin = group_min, xmax = group_max,
                ymin = ypos - 0.4, ymax = ypos + 0.4,
                fill = group_color)) +
  scale_fill_identity() +
  scale_y_reverse(breaks = modality_pos$ypos,
                     labels = modality_pos$modality,
                  expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  facet_grid(rows = vars(genotype_label)) +
  theme_classic(5) +
  theme(strip.text.y = element_text(angle = 0))

ggsave("genotype_group_barplots.pdf",
       group_bars,
       width = 4,
       height = 6
       )

# Counts for figure

txn_counts <- txn_rect_data %>%
  arrange(genotype_label) %>%
  select(genotype_label, genotype_n) %>%
  unique()

# Overlap stats for figure:
