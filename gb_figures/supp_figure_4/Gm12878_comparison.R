# This script contains platform comparisons of scATAC-seq methods using Gm12878 cells.
# 
# ## Data Sources
# 
# The datasets explored in this analysis are from the following sources:  
# 
# **bu**:  
# Buenrostro JD, Wu B, Litzenburger UM, Ruff D et al. Single-cell chromatin accessibility reveals principles of regulatory variation. Nature 2015 Jul 23;523(7561):486-90. PMID: 26083756
# 
# URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65360
# 
# **cu**:  
# Cusanovich DA, Daza R, Adey A, Pliner HA et al. Multiplex single cell profiling of chromatin accessibility by combinatorial cellular indexing. Science 2015 May 22;348(6237):910-4. PMID: 25953818
# 
# URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67446
# 
# **pl**:
# Pliner H, Packer J, McFaline-Figueroa et al. Cicero Predicts cis-Regulatory DNA Interactions from Single-Cell Chromatin Accessibility Data. Mol Cell 2018.
#
# URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109828
#
# **tx**:  
# 10X Genomics, 5k 1:1 mixture of fresh frozen human (GM12878) and mouse (A20) cells  
# 
# URL: https://support.10xgenomics.com/single-cell-atac/datasets/1.0.0/atac_v1_hgmm_5k
# 
# **gr**:  
# The present study (Graybuck, et al.).
# 
# Datasets bu, cu, and gr were aligned to mm10 using the same process described in Graybuck, et al. using raw data downloaded from GEO or generated using a Miseq (for gb).
# 
# Data from 10X genomics was obtained by directly downloading the fragments.tsv table from the 10X genomics website, followed by filtering to remove mouse cells and all barcodes with < 100 sequenced fragments.

## Read datasets

library(dplyr)
library(lowcat)
library(purrr)
library(GenomicAlignments)
library(ggplot2)
library(ggbeeswarm)
options(stringsAsFactors = F)

#Read the BAM files for the bu dataset:
bu_data_dir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE65360/raw/bam"

bu_data_files <- list.files(bu_data_dir, 
                            pattern = ".rmd.srt.bam$",
                            full.names = TRUE)

bu_fragments <- run_pe_to_frag_parallel(bu_data_files)

#Read the BAM files for the cu dataset:
cu_data_dir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE67446/raw/bam"

cu_data_files <- list.files(cu_data_dir, 
                            pattern = ".rmd.srt.bam$",
                            full.names = TRUE)

cu_fragments <- run_pe_to_frag_parallel(cu_data_files)


#Read the BAM files for the gr dataset:
gr_data_dir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/Gm12878_2017-03-29/raw/bam"

gr_data_files <- list.files(gr_data_dir,
                            pattern = ".rmd.srt.bam$",
                            full.names = TRUE)

gr_manifest <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/Gm12878_2017-03-29/human_sample_manifest.csv")

gr_manifest <- gr_manifest %>%
  filter(unique_fragments > 100)

gr_data_files <- gr_manifest$bam_file

gr_fragments <- run_pe_to_frag_parallel(gr_data_files)

# Read the demultiplexed regions for the pl dataset:
load("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE109828/SRR6652586_hg38_regions.rda")

pl_indexes <- read.delim("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE109828/GSM2970932_sciATAC_GM12878_HL60_indextable.txt.gz",
                         sep = "\t") %>%
  filter(putative_cell_type == "GM12878")

pl_fragments <- bam_regions[names(bam_regions) %in% pl_indexes$barcode]
rm(bam_regions)


#Read the pre-processed fragments for the tx dataset:
tx_data_dir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/Gm12878_10X_2018-11-15/"

load(file.path(tx_data_dir, "10X_bam_fragments.rda"))

tx_fragments <- bam_fragments
rm(bam_fragments)

## Alignment statistics and read stats

#Build sample tables for each dataset
bu_samples <- data.frame(file_name = bu_data_files) %>%
  mutate(sample_name = sub(".rmd.srt.bam","",file_name),
         sample_name = sub(".+/","",sample_name),
         dataset_short = "bu",
         dataset_full = "Buenrostro 2015")
names(bu_fragments) <- bu_samples$sample_name

cu_samples <- data.frame(file_name = cu_data_files) %>%
  mutate(sample_name = sub(".rmd.srt.bam","",file_name),
         sample_name = sub(".+/","",sample_name),
         dataset_short = "cu",
         dataset_full = "Cusanovich 2015")
names(cu_fragments) <- cu_samples$sample_name

gr_samples <- data.frame(file_name = gr_data_files) %>%
  mutate(sample_name = sub(".rmd.srt.bam","",file_name),
         sample_name = sub(".+/","",sample_name),
         dataset_short = "gr",
         dataset_full = "Graybuck 2019")
names(gr_fragments) <- gr_samples$sample_name

pl_samples <- data.frame(file_name = "SRR6652586_hg38.bam",
                         sample_name = names(pl_fragments),
                         dataset_short = "pl",
                         dataset_full = "Pliner 2018")

tx_samples <- data.frame(file_name = "atac_v1_hgmm_5k_fragments.tsv",
                         sample_name = names(tx_fragments),
                         dataset_short = "tx",
                         dataset_full = "10X Genomics 2018")

#Read alignment stats for the bu, cu, and gr datasets.
bu_stats <- read.table("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE65360/clean/stats/compiled_library_stats.tsv", header = T)
bu_samples <- left_join(bu_samples, bu_stats, by = c("sample_name" = "sample"))
bu_fragments <- bu_fragments[bu_samples$sample_name]

cu_stats <- read.table("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE67446/clean/stats/compiled_library_stats.tsv", header = T)
cu_samples <- left_join(cu_samples, cu_stats, by = c("sample_name" = "sample"))
cu_fragments <- cu_fragments[cu_samples$sample_name]

gr_stats <- read.table("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/Gm12878_2017-03-29//clean/stats/compiled_library_stats.tsv", header = T)
gr_samples <- left_join(gr_samples, gr_stats, by = c("sample_name" = "sample"))

# For pl, read the FASTQ-based barcode counts
pl_fastq_counts <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE109828/SRR6652586_1_fixed.counts_filtered1000.txt")
pl_mapped_counts <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/GSE109828/SRR6652586_hg38_aligned.counts_filtered1000.txt")

pl_unique <- map(pl_fragments, unique)

pl_stats <- data.frame(sample_name = pl_samples$sample_name,
                       total_reads = pl_fastq_counts$Freq[match(pl_samples$sample_name, pl_counts$bcs1)] * 2,
                       mapped_reads = pl_mapped_counts$Freq[match(pl_samples$sample_name, pl_mapped_counts$bcs1)],
                       unique_fragments = map_int(pl_unique, length)) %>%
  mutate(mapped_fragments = mapped_reads / 2,
         unique_reads = unique_fragments * 2,
         mapped_percent = round(mapped_reads/total_reads * 100, 2),
         unmapped_percent = 100 - round(mapped_reads/total_reads * 100, 2),
         unique_percent = round(unique_reads / mapped_reads * 100, 2),
         duplicate_percent = 100 - round(unique_reads / mapped_reads * 100, 2)) %>%
  select(sample_name, total_reads, mapped_reads, mapped_fragments, mapped_percent, unmapped_percent, unique_reads, unique_fragments, unique_percent, duplicate_percent)

pl_samples <- left_join(pl_samples, pl_stats, by = "sample_name")  

# Keep only unique fragments for downstream analysis
pl_fragments <- pl_unique
pl_fragments <- map(pl_fragments, sort)

#For tx, similar statistics are supplied by 10X:
raw_tx_stats <- read.csv("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/Gm12878_10X_2018-11-15/atac_v1_hgmm_5k_singlecell.csv", header = TRUE)
raw_tx_stats <- filter(raw_tx_stats, barcode %in% tx_samples$sample_name)

tx_stats <- data.frame(sample = raw_tx_stats$barcode,
                       total_reads = raw_tx_stats$total * 2,
                       mapped_reads = (raw_tx_stats$passed_filters + raw_tx_stats$duplicate)*2,
                       mapped_fragments = raw_tx_stats$passed_filters + raw_tx_stats$duplicate,
                       mapped_percent = round((raw_tx_stats$passed_filters + raw_tx_stats$duplicate) / (raw_tx_stats$total) * 100, 2),
                       unmapped_percent = 100 - round((raw_tx_stats$passed_filters + raw_tx_stats$duplicate) / (raw_tx_stats$total) * 100, 2),
                       unique_reads = (raw_tx_stats$passed_filters_hg19)*2,
                       unique_fragments = raw_tx_stats$passed_filters_hg19,
                       unique_percent = round(raw_tx_stats$passed_filters_hg19 / raw_tx_stats$total * 100, 2),
                       duplicate_percent = 100 - round(raw_tx_stats$passed_filters_hg19 / raw_tx_stats$total * 100, 2)
                       )

tx_samples <- left_join(tx_samples, tx_stats, by = c("sample_name" = "sample"))

# match the order of the fragments
tx_samples <- tx_samples[match(names(tx_fragments), tx_samples$sample_name),]


#Median fragment length may also be informative:

gr_med_width <- function(gr) {
  median(width(gr))
}

bu_samples$med_length <- map_dbl(bu_fragments, gr_med_width)
cu_samples$med_length <- map_dbl(cu_fragments, gr_med_width)
gr_samples$med_length <- map_dbl(gr_fragments, gr_med_width)
pl_samples$med_length <- map_dbl(pl_fragments, gr_med_width)
tx_samples$med_length <- map_dbl(tx_fragments, gr_med_width)

## ENCODE Overlaps

#Compute overlaps with ENCODE data:
#```{r Load ENCODE Peaks}
# ENCFF570AMI <- read.delim("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/controls/ENCFF570AMI.bed.gz", header = F, sep = "\t")
# 
# ENCFF570AMI_gr <- GRanges(seqnames = ENCFF570AMI[,1],
#                           IRanges(start = ENCFF570AMI[,2],
#                                   end = ENCFF570AMI[,3]))


if(!file.exists("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/controls/ENCFF773SPT.bed.gz")) {
  download.file("https://www.encodeproject.org/files/ENCFF773SPT/@@download/ENCFF773SPT.bed.gz",
                "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/controls/ENCFF773SPT.bed.gz")
}

ENCFF773SPT <- read.delim("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/controls/ENCFF773SPT.bed.gz", header = F, sep = "\t")

ENCFF773SPT_gr <- GRanges(seqnames = ENCFF773SPT[,1],
                          IRanges(start = ENCFF773SPT[,2],
                                  end = ENCFF773SPT[,3]))


#```

# 10X is aligned to hg19, not hg38. So, we'll need peaks for that.
if(!file.exists("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/controls/ENCFF206HYT.bed.gz")) {
  download.file("https://www.encodeproject.org/files/ENCFF773SPT/@@download/ENCFF206HYT.bed.gz",
                "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/controls/ENCFF206HYT.bed.gz")
}

ENCFF206HYT <- read.delim("//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/data/controls/ENCFF206HYT.bed.gz", header = F, sep = "\t")

ENCFF206HYT_gr <- GRanges(seqnames = ENCFF206HYT[,1],
                          IRanges(start = ENCFF206HYT[,2],
                                  end = ENCFF206HYT[,3]))

#Count overlaps for each dataset
bu_encode <- colSums(count_fragment_overlaps(bu_fragments,
                                             ENCFF773SPT_gr,
                                             aggregate = FALSE,
                                             binarize = FALSE))
bu_samples$encode_count <- bu_encode
bu_samples$encode_frac <- round(bu_encode / bu_samples$unique_fragments, 4)

cu_encode <- colSums(count_fragment_overlaps(cu_fragments,
                                             ENCFF773SPT_gr,
                                             aggregate = FALSE,
                                             binarize = FALSE))
cu_samples$encode_count <- cu_encode
cu_samples$encode_frac <- round(cu_encode / cu_samples$unique_fragments, 4)

gr_encode <- colSums(count_fragment_overlaps(gr_fragments,
                                             ENCFF773SPT_gr,
                                             aggregate = FALSE,
                                             binarize = FALSE))
gr_samples$encode_count <- gr_encode
gr_samples$encode_frac <- round(gr_encode / gr_samples$unique_fragments, 4)


pl_encode <- colSums(count_fragment_overlaps(pl_fragments,
                                             ENCFF773SPT_gr,
                                             aggregate = FALSE,
                                             binarize = FALSE))
pl_samples$encode_count <- pl_encode
pl_samples$encode_frac <- round(pl_encode / pl_samples$unique_fragments, 4)

tx_encode <- colSums(count_fragment_overlaps(tx_fragments,
                                             ENCFF206HYT_gr,
                                             aggregate = FALSE,
                                             binarize = FALSE))
tx_samples$encode_count <- tx_encode
tx_samples$encode_frac <- round(tx_encode / tx_samples$unique_fragments, 4)


## RefSeq TSS Overlaps

#Compute overlaps with RefSeq TSS Regions:

#tss <- get_tss_regions(expand = 5e3, genome = "hg19")
load("hg38_tss.rda")
hg38_tss <- hg38_tss %>%
  filter(!grepl("_",chr))
hg38_tss_gr <- bed_to_GRanges(hg38_tss)
hg38_tss_gr <- GenomicRanges::reduce(hg38_tss_gr)

load("hg19_tss.rda")
hg19_tss <- hg19_tss %>%
  filter(!grepl("_",chr))
hg19_tss_gr <- bed_to_GRanges(hg19_tss)
hg19_tss_gr <- GenomicRanges::reduce(hg19_tss_gr)


#Count overlaps for each dataset
bu_tss <- colSums(count_fragment_overlaps(bu_fragments,
                                          hg38_tss_gr,
                                          aggregate = FALSE,
                                          binarize = FALSE))
bu_samples$tss_count <- bu_tss
bu_samples$tss_frac <- round(bu_tss / bu_samples$unique_fragments, 4)

cu_tss <- colSums(count_fragment_overlaps(cu_fragments,
                                          hg38_tss_gr,
                                          aggregate = FALSE,
                                          binarize = FALSE))
cu_samples$tss_count <- cu_tss
cu_samples$tss_frac <- round(cu_tss / cu_samples$unique_fragments, 4)

gr_tss <- colSums(count_fragment_overlaps(gr_fragments,
                                          hg38_tss_gr,
                                          aggregate = FALSE,
                                          binarize = FALSE))
gr_samples$tss_count <- gr_tss
gr_samples$tss_frac <- round(gr_tss / gr_samples$unique_fragments, 4)

pl_tss <- colSums(count_fragment_overlaps(pl_fragments,
                                          hg38_tss_gr,
                                          aggregate = FALSE,
                                          binarize = FALSE))
pl_samples$tss_count <- pl_tss
pl_samples$tss_frac <- round(pl_tss / pl_samples$unique_fragments, 4)

tx_tss <- colSums(count_fragment_overlaps(tx_fragments,
                                          hg19_tss_gr,
                                          aggregate = FALSE,
                                          binarize = FALSE))
tx_samples$tss_count <- tx_tss
tx_samples$tss_frac <- round(tx_tss / tx_samples$unique_fragments, 4)


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

  if(roundall) {
    roundcols <- c("median","mean","max","sd","q25","q75","se","ci")
    datac[roundcols] <- round(datac[roundcols],3)
  }
  
  
  return(datac)
}


#First, we'll assemble all of the data into a single object:
all_samples <- rbind(bu_samples, cu_samples, gr_samples, pl_samples, tx_samples)

# Now, we can plot them together using the dataset_short column for grouping.
# 
# First, total fragments per sample:
total_reads_data <- select(all_samples, dataset_full,
                           sample_name, dataset_short, total_reads) %>%
  mutate(total_reads = log10(total_reads + 1))

total_reads_summary <- summarySE(data = total_reads_data,
                                 measurevar = "total_reads",
                                 groupvars = "dataset_short")

total_reads_plot <- ggplot() +
  geom_quasirandom(data = total_reads_data,
                   aes(x = dataset_short,
                       y = total_reads,
                       color = as.factor(dataset_full)),
                   size = 0.2) +
  geom_errorbar(data = total_reads_summary,
                aes(x = dataset_short,
                    ymin = q25,
                    ymax = q75)) +
  geom_point(data = total_reads_summary,
             aes(x = dataset_short,
                 y = median),
             color = "red",
             size = 1) +
  ylab("Log10(Total Reads + 1)") +
  xlab("Dataset") +
  ggtitle("Total Reads per Cell by Dataset") +
  theme_classic(5) +
  theme(legend.position = "none")

ggsave("total_reads_plot.pdf",
       total_reads_plot,
       width = 3, height = 2,
       useDingbats = FALSE)

#Unique Fragments per sample:
#```{r}
unique_fragments_data <- select(all_samples, dataset_full,
                           sample_name, dataset_short, unique_fragments) %>%
  mutate(unique_fragments = log10(unique_fragments + 1))

unique_fragments_summary <- summarySE(data = unique_fragments_data,
                                 measurevar = "unique_fragments",
                                 groupvars = "dataset_short")

unique_fragments_plot <- ggplot() +
  geom_quasirandom(data = unique_fragments_data,
                   aes(x = dataset_short,
                       y = unique_fragments,
                       color = as.factor(dataset_full)),
                   size = 0.2) +
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

#Percent Unique Fragments per sample:
unique_percent_data <- select(all_samples, dataset_full,
                           sample_name, dataset_short, unique_percent) %>%
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

# Median Fragment Length per sample:

med_length_data <- select(all_samples, dataset_full,
                           sample_name, dataset_short, med_length) %>%
  filter(!is.na(med_length))

med_length_summary <- summarySE(data = med_length_data,
                                 measurevar = "med_length",
                                 groupvars = "dataset_short")

median_length_plot <- ggplot() +
  geom_quasirandom(data = med_length_data,
                   aes(x = dataset_short,
                       y = med_length,
                       color = as.factor(dataset_full)),
                   size = 0.2) +
  geom_errorbar(data = med_length_summary,
                aes(x = dataset_short,
                    ymin = q25,
                    ymax = q75)) +
  geom_point(data = med_length_summary,
             aes(x = dataset_short,
                 y = median),
             color = "red",
             size = 1) +
  ylab("Median Fragment Length") +
  xlab("Dataset") +
  ggtitle("Median Fragment Length per Cell by Dataset") +
  theme_classic(5) +
  theme(legend.position = "none")

ggsave("median_length_plot.pdf",
       median_length_plot,
       width = 3, height = 2,
       useDingbats = FALSE)

#ENCODE Peak Overlap Fraction per sample:
encode_frac_data <- select(all_samples, dataset_full,
                           sample_name, dataset_short, encode_frac) %>%
  filter(!is.na(encode_frac))

encode_frac_summary <- summarySE(data = encode_frac_data,
                                 measurevar = "encode_frac",
                                 groupvars = "dataset_short")

encode_frac_plot <- ggplot() +
  geom_quasirandom(data = encode_frac_data,
                   aes(x = dataset_short,
                       y = encode_frac,
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

#Refseq TSS Overlap Fraction per sample:
tss_frac_data <- select(all_samples, dataset_full,
                           sample_name, dataset_short, tss_frac) %>%
  filter(!is.na(tss_frac))

tss_frac_summary <- summarySE(data = tss_frac_data,
                                 measurevar = "tss_frac",
                                 groupvars = "dataset_short")

tss_frac_plot <- ggplot() +
  geom_quasirandom(data = tss_frac_data,
                   aes(x = dataset_short,
                       y = tss_frac,
                       color = as.factor(dataset_full)),
                   size = 0.2) +
  geom_errorbar(data = tss_frac_summary,
                aes(x = dataset_short,
                    ymin = q25,
                    ymax = q75)) +
  geom_point(data = tss_frac_summary,
             aes(x = dataset_short,
                 y = median),
             color = "red",
             size = 1) +
  ylab("Refseq TSS Overlap Fraction") +
  xlab("Dataset") +
  ggtitle("Refseq TSS Overlap Fraction per Cell by Dataset") +
  theme_classic(5) +
  theme(legend.position = "none")

ggsave("tss_frac_plot.pdf",
       tss_frac_plot,
       width = 3, height = 2,
       useDingbats = FALSE)

#Plots of unique fragments vs ENCODE overlap:
all_samples$dataset_full <- as.factor(all_samples$dataset_full)

unique_vs_encode_plot <- ggplot() +
  geom_point(data = all_samples %>% filter(dataset_short == "tx"),
             aes(x = log10(unique_fragments + 1),
                 y = encode_frac,
                 color = dataset_full),
             size = 0.2) +
  geom_point(data = all_samples %>% filter(dataset_short != "tx"),
             aes(x = log10(unique_fragments + 1),
                 y = encode_frac,
                 color = dataset_full),
             size = 0.2) +
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

ggsave("unique_vs_encode_plot.pdf",
       unique_vs_encode_plot,
       width = 3, height = 3,
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

bu_frag_matrix <- fragment_size_matrix(bu_fragments)
cu_frag_matrix <- fragment_size_matrix(cu_fragments)
gr_frag_matrix <- fragment_size_matrix(gr_fragments)
pl_frag_matrix <- fragment_size_matrix(pl_fragments)
tx_frag_matrix <- fragment_size_matrix(tx_fragments)

bu_frag_df <- fragment_mat_to_df(bu_frag_matrix)
cu_frag_df <- fragment_mat_to_df(cu_frag_matrix)
gr_frag_df <- fragment_mat_to_df(gr_frag_matrix)
pl_frag_df <- fragment_mat_to_df(pl_frag_matrix)
tx_frag_df <- fragment_mat_to_df(tx_frag_matrix)

bu_frag_plot <- ggplot() +
  geom_ribbon(data = bu_frag_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "#A3A500",
              alpha = 0.4) +
  geom_line(data = bu_frag_df,
            aes(x = size,
                y = med),
            color = "#A3A500",
            size = 0.5) +
  scale_x_continuous("") +
  theme_classic(5)

cu_frag_plot <- ggplot() +
  geom_ribbon(data = cu_frag_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "#00BF7D",
              alpha = 0.4) +
  geom_line(data = cu_frag_df,
            aes(x = size,
                y = med),
            color = "#00BF7D",
            size = 0.5) +
  scale_x_continuous("") +
  theme_classic(5)

gr_frag_plot <- ggplot() +
  geom_ribbon(data = gr_frag_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "#00B0F6",
              alpha = 0.4) +
  geom_line(data = gr_frag_df,
            aes(x = size,
                y = med),
            color = "#00B0F6",
            size = 0.5) +
  scale_x_continuous("") +
  theme_classic(5)

pl_frag_plot <- ggplot() +
  geom_ribbon(data = pl_frag_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "#E76BF3",
              alpha = 0.4) +
  geom_line(data = pl_frag_df,
            aes(x = size,
                y = med),
            color = "#E76BF3",
            size = 0.5) +
  scale_x_continuous("") +
  theme_classic(5)

tx_frag_plot <- ggplot() +
  geom_ribbon(data = tx_frag_df,
              aes(x = size,
                  ymin = q25,
                  ymax = q75),
              fill = "#F8766D",
              alpha = 0.4) +
  geom_line(data = tx_frag_df,
            aes(x = size,
                y = med),
            color = "#F8766D",
            size = 0.5) +
  scale_x_continuous("") +
  theme_classic(5)

library(cowplot)

all_frag_plot <- plot_grid(bu_frag_plot,
                           cu_frag_plot,
                           gr_frag_plot,
                           pl_frag_plot,
                           tx_frag_plot,
                           ncol = 1)

ggsave("fragment_size_plots.pdf",
       all_frag_plot,
       width = 3, height = 4,
       useDingbats = FALSE)
