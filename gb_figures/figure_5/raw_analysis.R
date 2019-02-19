library(cocoframer)
library(reshape2)
library(dplyr)
vol_dims <- c(528,320,456)

ccf_arr <- read_aibs_nrrd("../../CCFv3_annotation_25.nrrd", dims = vol_dims)
ccf_anno <- read.table("../../CCFv3_annotation_ITKSNAP_labels_lg.txt")

mba_anno <- get_mba_ontology()
mba_anno <- flatten_mba_ontology(mba_anno)

names(ccf_anno) <- c("id","x","y","z","a","b","c","acronym")
ccf_anno <- ccf_anno %>%
  select(id, acronym)

ccf <- melt(ccf_arr)
names(ccf) <- c("x","y","z","id")

ccf <- ccf %>%
  left_join(ccf_anno)

ccf_summary <- ccf %>%
  group_by(id, acronym) %>%
  summarise(vox = n())

raw_image <- readBin("../../400225_Retroorbital_TissueCyte_mscRE4-FlpO/resampled_red.raw", "integer", size = 2, n = vol_dims[1] * vol_dims[2] * vol_dims[3])

raw_arr <- array(raw_image, dim = vol_dims)

dat <- melt(raw_arr)
names(dat) <- c("x","y","z","val")

dat$id <- ccf$id
dat$acronym <- ccf$acronym

dat <- dat %>%
  filter(!is.na(acronym))

cutoff <- 850

str_summary <- dat %>%
  group_by(id, acronym) %>%
  summarise(vox = n(),
            mean_val = mean(val),
            min_val = min(val),
            max_val = max(val),
            vox_gt_cut = sum(val > cutoff))

bg_sub_dat <- dat %>%
  mutate(val = val - cutoff)
bg_sub_dat$val[bg_sub_dat$val < 0] <- 0

bg_sub_str_summary <- bg_sub_dat %>%
  group_by(id, acronym) %>%
  summarise(vox = n(),
            sum_val = sum(val),
            mean_val = mean(val),
            min_val = min(val),
            max_val = max(val),
            mean_nz = mean(val[val > 0]),
            frac_nz = sum(val > 0)/n(),
            energy = sum(val) * frac_nz) %>%
  mutate(mean_nz = ifelse(is.na(mean_nz), 0, mean_nz),
         frac_nz = ifelse(is.na(frac_nz), 0, frac_nz))

bg_sub_str_summary <- bg_sub_str_summary %>%
  filter(acronym %in% mba_anno$acronym)

write.csv(bg_sub_str_summary,
          "background_subtracted_summary.csv",
          row.names = FALSE)
