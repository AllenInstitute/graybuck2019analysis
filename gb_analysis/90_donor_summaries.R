library(dplyr)
options(stringsAsFactors = F)

samples <- read.csv("../common/all_samples.csv")

# N samples
nrow(samples)

## Driver/reporter samples
driver_reporter_samples <- samples %>%
  filter(full_genotype != "Ai14(RCL-tdT)/wt") %>% # Retrograde
  filter(full_genotype != "wt/wt") # Retroorbital

# N driver/reporter samples
nrow(driver_reporter_samples)

# N driver/reporter genotypes
length(unique(driver_reporter_samples$full_genotype))

# N driver/reporter mice
length(unique(driver_reporter_samples$animal))

## Retrograde samples
retrograde_samples <- samples %>%
  filter(full_genotype == "Ai14(RCL-tdT)/wt")

# N retrograde
nrow(retrograde_samples)

# N retrograde locations
length(unique(retrograde_samples$inj_target))

# N retrograde mice
length(unique(retrograde_samples$animal))

## Retroorbital samples
retroorbital_samples <- samples %>%
  filter(full_genotype == "wt/wt")

# N retroorbital
nrow(retroorbital_samples)

# N retroorbital donors
length(unique(retroorbital_samples$animal))

## High quality samples
qc_samples <- samples %>%
  filter(frac_gt_250bp > 0.1,
         unique_fragments > 1e4,
         ENCODE_frac > 0.25)

# N QC Samples
nrow(qc_samples)

# QC UpSet
library(UpSetR)
samples_qc_binary <- samples %>%
  mutate(unique_fragments = as.numeric(unique_fragments <= 1e4),
         frac_gt_250bp = as.numeric(frac_gt_250bp <= 0.1),
         ENCODE_frac = as.numeric(ENCODE_frac <= 0.25)) %>%
  rowwise() %>%
  mutate(pass_qc = as.numeric(unique_fragments + frac_gt_250bp + ENCODE_frac == 0)) %>%
  ungroup() %>%
  select(sample_id, unique_fragments, frac_gt_250bp, ENCODE_frac, pass_qc) %>%
  as.data.frame()

pdf("qc_upset.pdf",
    width = 3.5, height = 2,
    useDingbats = FALSE)
upset(samples_qc_binary, 
      sets = c("unique_fragments","frac_gt_250bp","ENCODE_frac","pass_qc"),
      text.scale = 0.5)
dev.off()
