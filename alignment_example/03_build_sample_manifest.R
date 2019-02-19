library(dplyr)
library(plater)

donors <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/donors.csv") %>%
  select(-notes,-n_cells)

data_dir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/2018-11-08/"

bam_files <- list.files(paste0(data_dir,"raw/bam/"), pattern = "rmd.srt.bam$")
bam_stats <- read.table(paste0(data_dir,"clean/stats/compiled_library_stats.tsv"), header = TRUE)

plate_df <- read_plate(file.path(data_dir,"sample.info.plater_batch41.csv")) %>%
  mutate(sample_id = paste0("20181108_",Wells))

keep_columns <- c("sample_id","i5_index","i7_index","animal")

sample_manifest <- plate_df[,keep_columns] %>%
  left_join(donors) %>%
  rowwise() %>%
  mutate(bam_file = paste0(data_dir,"raw/bam/",bam_files[grepl(paste0("N",i7_index,"-","N",i5_index), bam_files)])) %>%
  mutate(sample = sub(".+/","",sub("\\..+","",bam_file))) %>%
  left_join(bam_stats) %>%
  select(-sample) %>%
  filter(!is.na(i5_index))

write.csv(sample_manifest,
          file.path(data_dir,"sample_manifest.csv"),
          quote = F, row.names = F)

