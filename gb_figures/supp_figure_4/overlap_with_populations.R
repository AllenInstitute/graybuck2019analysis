library(dplyr)
library(purrr)
library(lowcat)

adir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/scatac_analysis/doc/paper/2018-02-15_paper_dataset/"
popdir <- "//allen/programs/celltypes/workgroups/mct-t200/T502/atac_analysis/macs2_3.2M/"

# Load GRanges for all bams
load(file.path(adir,"bam_fragments.rda"))
# Load QC-passed samples
samples <- read.csv(file.path(adir,"gt1e4_samples.csv"), row.names = 1)
# Filter bams
bam_fragments <- bam_fragments[samples$sample_id]

# Load population MACS2 peaks
peak_files <- list.files(popdir, pattern = ".narrowPeak", full.names = T)
peak_files <- peak_files[!grepl("genomic",peak_files)]

peak_names <- c("cux2_rep1","cux2_rep2","cux2_rep3",
                "gad2_rep1","gad2_rep2","gad2_rep3",
                "ntsr1_rep1","ntsr1_rep2","ntsr1_rep3",
                "rbp4_rep1","rbp4_rep2",
                "scnn1a_rep1","scnn1a_rep2","scnn1a_rep3",
                "mes_rep1","mes_rep2","mes_rep3",
                "camk2a_rep1","camk2a_rep2",
                "pvalb_rep1","pvalb_rep2",
                "vip_rep1","vip_rep2",
                "rbp4_rep3")

peak_beds <- map(peak_files, 
                 function(x) {
                   bed <- read.table(x)
                   names(bed) <- c("chr","start","end","name","score","strand","fc","nlog10p","nlog10q","rel_summit_pos")
                   bed$strand <- "+"
                   bed
                 } )
names(peak_beds) <- peak_names

peak_GR <- map(peak_beds, bed_to_GRanges)
names(peak_GR) <- peak_names

calculate_GR_counts <- function(peak_list, bam_list) {
  
  counts <- matrix(nrow=length(bam_list),ncol=length(peak_list))
  
  for(b in 1:length(bam_list)) {
    print(paste0("Running ",b," of ",length(bam_list)))
    
    bam <- bam_list[[b]]
    
    for(p in 1:length(peak_list)) {
      counts[b,p] <- sum(countOverlaps(bam,peak_list[[p]]))
      
    }
  }
  
  counts <- as.data.frame(counts)
  names(counts) <- names(peak_list)
  row.names(counts) <- names(bam_list)
  
  return(counts)
  
}

peak_counts <- calculate_GR_counts(peak_GR,
                                   bam_fragments)

n_reads <- map_int(bam_fragments, length)

pc_mat <- as.matrix(peak_counts)
pf_mat <- apply(pc_mat, 2, function(x) x/n_reads)

heatmap(pf_mat)

pn_mat <- apply(pf_mat, 2, function(x) {
  if(max(x) > 0) { x/max(x) } else {
    rep(0, length(x))
    
  }
  }
  )
heatmap(pn_mat)

# Jaccard distances based on overlaps

pj_mat <- pc_mat
peak_n <- map_int(peak_GR, length)
bam_n <- map_int(bam_fragments, length)

for(p in 1:ncol(pj_mat)) {
  for(b in 1:nrow(pj_mat)) {
    pj_mat[b,p] <- pj_mat[b,p]/(bam_n[b] + peak_n[p])
  }
}

heatmap(pj_mat,
        hclustfun = function(x) hclust(x, method = "ward.D2"),
        Colv = NA,
        Rowv = NA)

### Interesting, but not useful.
library(Rtsne)

peak_tsne <- Rtsne(pj_mat[,!grepl("mes|pvalb|vip|camk2a", colnames(pj_mat))], check_duplicates = FALSE, perplexity = 10)

tsne_df <- data.frame(x = peak_tsne$Y[,1],
                      y = peak_tsne$Y[,2])

ggplot() +
  geom_point(data = tsne_df,
             aes(x = x,
                 y = y,
                 color = samples$source_label))
