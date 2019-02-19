# ## Data Sources
# 
# The datasets explored in this analysis are from the following sources:  
# 
# **cu**:  
# Cusanovich DA, Hill AJ, Aghamirzaie D, Daza RM et al. 
# A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. 
# Cell 2018 Aug 23;174(5):1309-1324.e18. PMID: 30078704
# 
# URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111586
# 
# Specifically, the sample from Prefrontal cortex provided by GSM3034633:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3034633
#
# FASTQs were obtained directly from Andrew Hill and Darren Cusanovich because of barcode removal
# by the files in the Sequence Read Archive.
#
# **pr**:  
# Preissl S, Fang R, Huang H, Zhao Y et al. 
# Single-nucleus analysis of accessible chromatin in developing mouse forebrain reveals 
# cell-type-specific transcriptional regulation. 
# Nat Neurosci 2018 Mar;21(3):432-439. PMID: 29434377
# 
# URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100033
#
# Specifically, the sample for P65 Forebrain scATAC-seq provided by GSM2668124:
# https://www.ncbi.nlm.nih.gov/sra/?term=SRR6768122
# 
# **gr**:  
# The present study (Graybuck, et al.).
# 
# Datasets cu, pr, and gr were aligned to mm10 using the same process described in Graybuck, et al. 
# 

## Read datasets
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggExtra)
library(lowcat)
library(purrr)
library(GenomicAlignments)
options(stringsAsFactors = F)

#Read the stats for cu
cu_samples <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2019-01-11_GSE111586_analysis/all_samples.csv")
cu_samples <- cu_samples %>%
  mutate(dataset_full = "Cusanovich_2018",
         dataset_short = "cu")

#Read the stats for the gr dataset:
gr_samples <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/all_samples.csv")
gr_samples <- gr_samples %>%
  filter(miseq_batch < 20181103) %>%
  mutate(dataset_full = "Graybuck_2019",
         dataset_short = "gr")

# Read stats for the pr dataset:
pr_samples <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2019-01-10_GSE100033_analysis/all_samples.csv")
pr_samples <- pr_samples %>%
  mutate(dataset_full = "Preissl_2018",
         dataset_short = "pr")

## Plots of alignment and overlap statistics

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE, roundall = F) {
  require(dplyr)
  # This does the summary. For each group return a vector with
  # N, mean, and sd
  
  names(data)[names(data) == measurevar] <- "measurevar"
  
  datac <- data %>%
    select(one_of(groupvars,"measurevar")) %>%
    filter(ifelse(na.rm == T, !is.na(measurevar), T)) %>%
    mutate(measurevar = as.numeric(measurevar)) %>%
    group_by_(c(groupvars)) %>%
    summarise(N = n(),
              median = median(measurevar),
              mean = mean(measurevar),
              max = max(measurevar),
              sd = ifelse(N == 1, 0, sd(measurevar)),
              q25 = as.numeric(quantile(measurevar, 0.25)),
              q75 = as.numeric(quantile(measurevar, 0.75))) %>%
    mutate(se = sd/sqrt(N))
  #%>%
  #  mutate(ci =  se * qt(conf.interval/2 + 0.5, N-1))
  
  
  if(roundall) {
    roundcols <- c("median","mean","max","sd","q25","q75","se","ci")
    datac[roundcols] <- round(datac[roundcols],3)
  }
  
  # datac <- datac %>%
  #   mutate(xpos = 1:n())
  
  return(datac)
}

#First, we'll assemble all of the data into a single object:
common_cols <- intersect(names(cu_samples), names(gr_samples))
common_cols <- intersect(common_cols, names(pr_samples))

all_samples <- rbind(cu_samples[, common_cols], 
                     gr_samples[, common_cols], 
                     pr_samples[, common_cols])

all_samples <- all_samples %>%
  mutate(unique_percent = unique_fragments / mapped_fragments * 100)
#```

# Now, we can plot them together using the dataset_short column for grouping.
# 
# First, total fragments per sample:
#```{r}
mapped_fragments_data <- select(all_samples, dataset_full,
                           sample_id, dataset_short, mapped_fragments) %>%
  mutate(mapped_fragments = log10(mapped_fragments + 1))

mapped_fragments_summary <- summarySE(data = mapped_fragments_data,
                                 measurevar = "mapped_fragments",
                                 groupvars = "dataset_short")

mapped_fragments_plot <- ggplot() +
  geom_quasirandom(data = mapped_fragments_data,
                   aes(x = dataset_short,
                       y = mapped_fragments,
                       color = as.factor(dataset_full)),
                   size = 0.1) +
  geom_errorbar(data = mapped_fragments_summary,
                aes(x = dataset_short,
                    ymin = q25,
                    ymax = q75)) +
  geom_point(data = mapped_fragments_summary,
             aes(x = dataset_short,
                 y = median),
             color = "red",
             size = 1) +
  ylab("Log10(Mapped Fragments + 1)") +
  xlab("Dataset") +
  ggtitle("Mapped Fragments per Cell by Dataset") +
  theme_classic(5) +
  theme(legend.position = "none")

ggsave("mapped_fragments_plot.pdf",
       mapped_fragments_plot,
       width = 3, height = 2,
       useDingbats = FALSE)
#```

#Unique Fragments per sample:
#```{r}
unique_fragments_data <- select(all_samples, dataset_full,
                           sample_id, dataset_short, unique_fragments) %>%
  mutate(unique_fragments = log10(unique_fragments + 1))

unique_fragments_summary <- summarySE(data = unique_fragments_data,
                                 measurevar = "unique_fragments",
                                 groupvars = "dataset_short")

unique_fragments_plot <- ggplot() +
  geom_quasirandom(data = unique_fragments_data,
                   aes(x = dataset_short,
                       y = unique_fragments,
                       color = as.factor(dataset_full)),
                   size = 0.1) +
  geom_errorbar(data = unique_fragments_summary,
                aes(x = dataset_short,
                    ymin = q25,
                    ymax = q75)) +
  geom_point(data = unique_fragments_summary,
             aes(x = dataset_short,
                 y = median),
             color = "red",
             size = 1) +
  ylab("Log10(Unique Fragments + 1)") +
  xlab("Dataset") +
  ggtitle("Unique Fragments per Cell by Dataset") +
  theme_classic(5) +
  theme(legend.position = "none")

ggsave("unique_fragments_plot.pdf",
       unique_fragments_plot,
       width = 3, height = 2,
       useDingbats = FALSE)
#```

#Percent Unique Fragments per sample:
#```{r}
unique_percent_data <- select(all_samples, dataset_full,
                           sample_id, dataset_short, unique_percent) %>%
  filter(!is.na(unique_percent))

unique_percent_summary <- summarySE(data = unique_percent_data,
                                 measurevar = "unique_percent",
                                 groupvars = "dataset_short")

unique_percent_plot <- ggplot() +
  geom_quasirandom(data = unique_percent_data,
                   aes(x = dataset_short,
                       y = unique_percent,
                       color = as.factor(dataset_full)),
                   size = 0.2) +
  geom_errorbar(data = unique_percent_summary,
                aes(x = dataset_short,
                    ymin = q25,
                    ymax = q75)) +
  geom_point(data = unique_percent_summary,
             aes(x = dataset_short,
                 y = median),
             color = "red",
             size = 1) +
  ylab("Percent Unique Fragments") +
  xlab("Dataset") +
  ggtitle("Percent Unique Fragments per Cell by Dataset") +
  theme_classic(5) +
  theme(legend.position = "none")

ggsave("unique_percent_plot.pdf",
       unique_percent_plot,
       width = 3, height = 2,
       useDingbats = FALSE)
#```


#ENCODE Peak Overlap Fraction per sample:
#```{r}
encode_frac_data <- select(all_samples, dataset_full,
                           sample_id, dataset_short, ENCODE_frac) %>%
  filter(!is.na(ENCODE_frac))

encode_frac_summary <- summarySE(data = encode_frac_data,
                                 measurevar = "ENCODE_frac",
                                 groupvars = "dataset_short")

encode_frac_plot <- ggplot() +
  geom_quasirandom(data = encode_frac_data,
                   aes(x = dataset_short,
                       y = ENCODE_frac,
                       color = as.factor(dataset_full)),
                   size = 0.2) +
  geom_errorbar(data = encode_frac_summary,
                aes(x = dataset_short,
                    ymin = q25,
                    ymax = q75)) +
  geom_point(data = encode_frac_summary,
             aes(x = dataset_short,
                 y = median),
             color = "red",
             size = 1) +
  ylab("ENCODE Peak Overlap Fraction") +
  xlab("Dataset") +
  ggtitle("ENCODE Peak Overlap Fraction per Cell by Dataset") +
  theme_classic(5) +
  theme(legend.position = "none")

ggsave("encode_frac_plot.pdf",
       encode_frac_plot,
       width = 3, height = 2,
       useDingbats = FALSE)
#```



#Plots of unique fragments vs ENCODE overlap:
#```{r unique vs ENCODE plot}
all_samples$dataset_full <- as.factor(all_samples$dataset_full)

unique_vs_encode_plot <- ggplot() +
  geom_point(data = all_samples,
             aes(x = log10(unique_fragments + 1),
                 y = ENCODE_frac,
                 color = dataset_full),
             size = 0.1) +
  geom_hline(aes(yintercept = 0.25),
             size = 0.1,
             linetype = "dashed") +
  geom_vline(aes(xintercept = 4),
             size = 0.1,
             linetype = "dashed") +
  theme_classic(5)  +
  scale_x_continuous("log10(Unique Fragments + 1)") +
  scale_y_continuous("ENCODE Peak Overlap Fraction") +
  theme(legend.position = "none")

unique_vs_encode_mp <- ggMarginal(unique_vs_encode_plot,
                                  groupColour = TRUE,
                                  groupFill = TRUE)

ggsave("unique_vs_encode_plot.pdf",
       unique_vs_encode_mp,
       width = 4, height = 4,
       useDingbats = FALSE)
#```

## Fragment Size Analysis

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

fragment_mat_to_df <- function(fragment_mat) {
  data.frame(size = as.numeric(rownames(fragment_mat)),
             med = rowMedians(fragment_mat),
             q25 = apply(fragment_mat, 1, function(x) {
               quantile(x)[2]
             }),
             q75 = apply(fragment_mat, 1, function(x) {
               quantile(x)[4]
             }))
}

load("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE111586/SC_atac.PreFrontalCortex_62216_bam_regions.rda")
bam_regions <- bam_regions[cu_samples$sample_id]
bam_regions <- map(bam_regions, unique)
cu_frag_matrix <- fragment_size_matrix(bam_regions)

load("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/common/bam_fragments.rda")
gr_frag_matrix <- fragment_size_matrix(bam_fragments)

load("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE100033/SRR6768122_bam_regions.rda")
bam_regions <- map(bam_regions, unique)
bam_regions_lens <- map_int(bam_regions, length)
bam_regions <- bam_regions[bam_regions_lens > 1000]
names(bam_regions) <- pr_samples$sample_id
pr_frag_matrix <- fragment_size_matrix(bam_regions)

cu_frag_df <- fragment_mat_to_df(cu_frag_matrix)
gr_frag_df <- fragment_mat_to_df(gr_frag_matrix)
pr_frag_df <- fragment_mat_to_df(pr_frag_matrix)


cu_frag_plot <- ggplot() +
  geom_ribbon(data = cu_frag_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "#00BFC4",
              alpha = 0.4) +
  geom_line(data = cu_frag_df,
            aes(x = size,
                y = med),
            color = "#00BFC4",
            size = 0.5) +
  scale_x_continuous("") +
  theme_classic(5)

gr_frag_plot <- ggplot() +
  geom_ribbon(data = gr_frag_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "#C77CFF",
              alpha = 0.4) +
  geom_line(data = gr_frag_df,
            aes(x = size,
                y = med),
            color = "#C77CFF",
            size = 0.5) +
  scale_x_continuous("") +
  theme_classic(5)

pr_frag_plot <- ggplot() +
  geom_ribbon(data = pr_frag_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "#F8766D",
              alpha = 0.4) +
  geom_line(data = pr_frag_df,
            aes(x = size,
                y = med),
            color = "#F8766D",
            size = 0.5) +
  scale_x_continuous("") +
  theme_classic(5)

library(cowplot)

all_frag_plot <- plot_grid(cu_frag_plot,
                           gr_frag_plot,
                           pr_frag_plot,
                           ncol = 1)

ggsave("fragment_size_plots.pdf",
       all_frag_plot,
       width = 3, height = 3,
       useDingbats = FALSE)
