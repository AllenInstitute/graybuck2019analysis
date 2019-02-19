devtools::install_github("AllenInstitute/scrattch.hicat", ref = "dev")

library(scrattch.hicat)
library(scrattch.io)
library(dplyr)
library(purrr)
library(Matrix)
options(stringsAsFactors = F)

load("../common/tss_regions.rda")

tss_regions$start <- tss_regions$start + 2e4
tss_regions$end <- tss_regions$end - 2e4
tss_regions <- tss_regions[,-5]

tome <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/tomes/facs/mouse_V1_ALM_20180520/transcrip.tome"

data <- read_tome_dgCMatrix(tome,"/data/exon")
data <- Matrix::t(data)
sample_names <- read_tome_sample_names(tome)
gene_names <- read_tome_gene_names(tome)

colnames(data) <- sample_names
rownames(data) <- gene_names

anno <- read_tome_anno(tome)

anno <- anno %>%
  filter(cluster_id < 134) %>%
  filter(!grepl("ALM",cluster_label))

data <- data[,anno$sample_name]

norm_data <- log2(data + 1)

# group-level comparisons
grouping <- read.csv("../common/cluster_grouping_for_tracks.csv")
names(grouping) <- sub("pred_","",names(grouping))

anno <- left_join(anno, grouping)

group_anno <- anno %>%
  select(group_id, group_label) %>%
  unique() %>%
  arrange(group_id)

group_cl <- as.factor(anno$group_id)
names(group_cl) <- anno$sample_name

group_markers <- map(group_anno$group_id,
                     function(x) {
                       print(paste("group",x))
                       group_specific_markers(x,
                                              norm_data,
                                              group_cl,
                                              de.param = de_param(),
                                              n.markers = 20)
                     })
names(group_markers) <- group_anno$group_label

group_markers_df <- map_dfr(1:length(group_markers),
                            function(x) {
                              df <- group_markers[[x]]
                              df$group_label <- names(group_markers)[x]
                              names(df)[1] <- "gene_name"
                              df %>% 
                                left_join(tss_regions, 
                                          by = c("gene_name" = "name")) %>% 
                                select(group_label,
                                       gene_name, chr, start, end, strand,
                                       everything())
                            })

write.csv(group_markers_df, "../common/rnaseq_group_markers.csv", row.names = F)

# subclass-level
subclass_anno <- anno %>%
  select(subclass_id, subclass_label) %>%
  unique() %>%
  arrange(subclass_id)

subclass_cl <- as.factor(anno$subclass_id)
names(subclass_cl) <- anno$sample_name

possible_markers <- possibly(group_specific_markers, NULL)

subclass_markers <- map(subclass_anno$subclass_id,
                        function(x) {
                          print(paste("subclass",x))
                          possible_markers(x,
                                                 norm_data,
                                                 subclass_cl,
                                                 de.param = de_param(),
                                                 n.markers = 20)
                        })
names(subclass_markers) <- subclass_anno$subclass_label

subclass_markers_df <- map_dfr(1:length(subclass_markers),
                            function(x) {
                              df <- subclass_markers[[x]]
                              if(!is.null(df)) {
                                df$subclass_label <- names(subclass_markers)[x]
                                names(df)[1] <- "gene_name"
                                df %>% 
                                  left_join(tss_regions, 
                                            by = c("gene_name" = "name")) %>% 
                                  select(subclass_label,
                                         gene_name, chr, start, end, strand,
                                         everything())
                              }
                              
                            })

write.csv(subclass_markers_df, "../common/rnaseq_subclass_markers.csv", row.names = F)

# cluster-level
cluster_anno <- anno %>%
  select(cluster_id, cluster_label) %>%
  unique() %>%
  arrange(cluster_id)

cluster_cl <- as.factor(anno$cluster_id)
names(cluster_cl) <- anno$sample_name

cluster_markers <- map(cluster_anno$cluster_id,
                        function(x) {
                          print(paste("cluster",x))
                          possible_markers(x,
                                                 norm_data,
                                                 cluster_cl,
                                                 de.param = de_param(),
                                                 n.markers = 20)
                        })
names(cluster_markers) <- cluster_anno$cluster_label

cluster_markers_df <- map_dfr(1:length(cluster_markers),
                            function(x) {
                              df <- cluster_markers[[x]]
                              if(!is.null(df)) {
                                df$cluster_label <- names(cluster_markers)[x]
                                names(df)[1] <- "gene_name"
                                df %>% 
                                  left_join(tss_regions, 
                                            by = c("gene_name" = "name")) %>% 
                                  select(cluster_label,
                                         gene_name, chr, start, end, strand,
                                         everything())
                              }
                             
                            })

write.csv(cluster_markers_df, "../common/rnaseq_cluster_markers.csv", row.names = F)
