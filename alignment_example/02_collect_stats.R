library("dplyr")
library("ggplot2")
options(stringsAsFactors=F)

dataset <- "2018-11-08"

# Import all configurations from the pipeline
# getwd()
# system("source ./00_config.sh")
# 
# stats_dir <- Sys.getenv("CONFIG_STATS_DIRECTORY")
# inserts_dir <- Sys.getenv("CONFIG_INSERTS_DIRECTORY")
# stats_results_dir <- Sys.getenv("CONFIG_STATS_RES_DIR")

stats_dir <- paste0("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/",dataset,"/clean/stats/")
inserts_dir <- paste0("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/",dataset,"/clean/inserts/")
stats_results_dir <- paste0("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/",dataset,"/clean/stats/")


stats_files <- list.files(stats_dir, pattern = "stats")
inserts_files <- list.files(inserts_dir, pattern = ".txt")

compiled <- data.frame(sample = character(),
                       total_reads = numeric(),
                       mapped_reads = numeric(),
                       mapped_fragments = numeric(),
                       mapped_percent = numeric(),
                       unmapped_percent = numeric(),
                       unique_reads = numeric(),
                       unique_fragments = numeric(),
                       unique_percent = numeric(),
                       duplicate_percent = numeric())

for(stats_file in stats_files) {
  
  stats <- read.delim(paste0(stats_dir,stats_file),header=F,sep=" ")
  stats[,1] <- as.numeric(stats[,1])
  
  results <- data.frame(sample = sub(".stats","",stats_file),
                        total_reads = stats[1,1],
                        mapped_reads = stats[3,1],
                        mapped_fragments = stats[3,1] / 2,
                        mapped_percent = round( (stats[3,1] / stats[1,1]) * 100, 2),
                        unmapped_percent = 100 - round( (stats[3,1] / stats[1,1]) * 100, 2),
                        unique_reads = stats[13,1],
                        unique_fragments = stats[17,1],
                        unique_percent = round( stats[13,1] / stats[3,1] * 100, 2),
                        duplicate_percent = 100 - round( stats[13,1] / stats[3,1] * 100, 2)
  )
  
  compiled <- rbind(compiled,results)
  
}

mapping_plot_stats <- compiled %>%
  mutate(unmapped_ypos = (total_reads + mapped_reads)/2,
         duplicate_ypos = (mapped_reads + unique_reads)/2,
         unique_ypos = unique_reads/2)

mapping_colors <- c("Unmapped" = "gray50",
                    "Duplicate" = "darkblue",
                    "Unique" = "dodgerblue")

mapping_plot <- ggplot() +
  geom_bar(data = compiled,
           aes(x = sample,
               y = total_reads,
               fill = "Unmapped"),
           stat = "identity") +
  geom_bar(data = compiled,
           aes(x = sample,
               y = mapped_reads,
               fill = "Duplicate"),
           stat = "identity") +
  geom_bar(data = compiled,
           aes(x = sample,
               y = unique_reads,
               fill = "Unique"),
           stat = "identity") +
  scale_fill_manual(name = "Category", values = mapping_colors) +
  geom_text(data = mapping_plot_stats,
            aes(x = sample,
                y = unmapped_ypos,
                label = paste0(unmapped_percent,"%")),
            color = "white",
            size = 2) +
  geom_text(data = mapping_plot_stats,
            aes(x = sample,
                y = duplicate_ypos,
                label = paste0(duplicate_percent,"%")),
            color = "white",
            size = 2) +
  geom_text(data = mapping_plot_stats,
            aes(x = sample,
                y = unique_ypos,
                label = paste0(unique_percent,"%")),
            color = "white",
            size = 2) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))

ggsave("mapping_summary_plot.pdf", mapping_plot, path= stats_results_dir, width = 17.5, height = 16)

write.table(compiled,paste(stats_dir,"compiled_library_stats.tsv", sep="/"),sep="\t",quote=F,row.names=F)
