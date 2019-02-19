library(dplyr)
library(Matrix)
library(matrixStats)
library(purrr)

library(scrattch.hicat)
library(scrattch.io)

options(stringsAsFactors = FALSE)


# Load Data

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/R_Object/SparseMatrix/20190204_RSC-184-192_mouse_star2.0_exon_sparse.Rdata")
new_exon <- exon
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/R_Object/SparseMatrix/20181119_RSC-004-183_mouse_star2.0_exon_sparse.Rdata")
exon <- cbind(exon, new_exon)


# Load Sample annotations

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/R_Object/SparseMatrix/20190204_RSC-184-192_mouse_star2.0_samp.dat.Rdata")
new_samp.dat <- samp.dat
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/R_Object/SparseMatrix/20181119_RSC-004-183_mouse_star2.0_samp.dat.Rdata")
samp.dat[sapply(samp.dat, is.factor)] <- lapply(samp.dat[sapply(samp.dat, is.factor)], 
                                                as.character)

samp.dat <- rbind(samp.dat, new_samp.dat)


# Select retroorbital donors:  
# mscRE1-SYFP2 L5: c("362845","362848")  
# mscRE4-SYFP2 L5: c("362849","362852")  
# mscRE4-FlpO;Ai65F L5: "389708" 
# mscRE4-FlpO;Ai65F All: "422331"
# mscRE10-FlpO;Ai65F All: "422566"
# mscRE16-FlpO;Ai65F All: "422569"
# 
# Select stereotaxic donors: "413054","413089","412835","413247"
# mscRE4-EGFP: "413054","413089"
# mscRE16-EGFP: "412835","413247"


viruses <- data.frame(virus = c("mscRE1-SYFP2","mscRE1-SYFP2",
                                "mscRE4-SYFP2","mscRE4-SYFP2",
                                "mscRE4-FlpO", "mscRE4-FlpO",
                                "mscRE10-FlpO", "mscRE16-FlpO",
                                "mscRE4-EGFP","mscRE4-EGFP",
                                "mscRE16-EGFP","mscRE16-EGFP"),
                      injection = c(rep("ro",8),
                                    rep("st",4)),
                      dissection = c("L5","L5",
                                     "L5","L5",
                                     "L5","All",
                                     "All","All",
                                     "All","All",
                                     "All","All"),
                      external_donor_name = c("362845","362848",
                                              "362849","362852",
                                              "389708","422331",
                                              "422566","422569",
                                              "413054","413089",
                                              "412835","413247"))

ro_samples <- samp.dat %>%
  filter(external_donor_name %in% viruses$external_donor_name) %>%
  left_join(viruses)

table(ro_samples$virus)

# Experiments - each hemisphere in stereotaxic experiments was injected with
# a different volume, so should be treated separately.

ro_samples <- ro_samples %>%
  mutate(experiment = ifelse(injection == "st",
                             paste(virus, injection, hemisphere_name, sep = "_"),
                             paste(virus, injection, sep = "_")))
table(ro_samples$experiment, ro_samples$external_donor_name)

#Select data for RO samples

ro_exon <- exon[,ro_samples$exp_component_name]

#Convert to log2(CPM + 1) for mapping
ro_data <- log2(cpm(ro_exon) + 1)
ro_data <- as.matrix(ro_data)

# Positive control - Anterograde donors into Npr3-IRES2-Cre VISp:  
#   AAV-FLEX-EGFP: "287122"
ant_samples <- samp.dat %>%
  filter(external_donor_name == "287122")

#Select and convert anterograde data
ant_exon <- exon[,ant_samples$exp_component_name]
ant_data <- log2(cpm(ant_exon) + 1)
ant_data <- as.matrix(ant_data)

#Unload huge exon matrices
rm(exon)
rm(new_exon)


#Load VISp reference data for mapping:
pdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/V1_ALM/process_24411/"
load(file.path(pdir,"cl.clean.rda"))
load(file.path(pdir,"cl.df.rda"))
load(file.path(pdir,"norm.dat.rda"))
load(file.path(pdir,"samp.dat.rda"))
load(file.path(pdir,"cl.med.rda"))
load(file.path(pdir,"V1.cl.rda"))

#Filter RO data to match norm.dat
common_genes <- intersect(rownames(ro_data), rownames(norm.dat))
ro_data <- ro_data[common_genes,]
norm.dat <- norm.dat[common_genes,]


#Load immune-related gene set
immune.genes <- scan("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/V1_ALM/process_21729/WGCNA_result_anterograde/interferon_type_I",what="character")



#Filter genes and perform the mapping using map_sampling() from scrattch.hicat
select.markers <- intersect(V1.markers, row.names(ro_data))
select.markers <- select.markers[rowMaxs(ro_data[select.markers,]) > 0]

map.result <- map_sampling(norm.dat[,V1.cells], 
                           droplevels(V1.cl[V1.cells]), 
                           ro_data, 
                           markers = select.markers)


#Save mapping results

save(map.result, file = "VISp_mapping_results.rda")


#Anterograde mapping

ant_data <- ant_data[common_genes,]
ant.map.result <- map_sampling(norm.dat[,V1.cells], 
                               droplevels(V1.cl[V1.cells]), 
                               ant_data, 
                               markers = select.markers)

save(ant.map.result, file = "Anterograde_VISp_mapping_results.rda")


# For comparisons, we'll remove the retrograde cells
V1.samp.dat <- samp.dat[match(V1.cells, samp.dat$exp_component_name),]

V1.samp.dat <- V1.samp.dat %>%
  filter(Injection_type != "retrograde")
V1.keep.cl <- V1.cl[names(V1.cl) %in% V1.samp.dat$exp_component_name]


#Run the analysis for each experiment category.

map.df <- map.result$map.df
experiments <- unique(ro_samples$experiment)
V1_table <- table(V1.keep.cl)
V1_keep <- V1_table[V1_table >= 10]

source("helper_functions.R")
library(ggplot2)
library(ggrepel)

de_results <- map(experiments,
                  function(experiment) {
                    
                    exp_samples <- ro_samples[ro_samples$experiment == experiment,]
                    exp_map.df <- map.df[exp_samples$exp_component_name,]
                    
                    map_table <- table(exp_map.df$pred.cl)
                    map_cl <- map_table[map_table >= 10]
                    
                    keep_cl <- intersect(names(map_cl), names(V1_keep))
                    
                    fg_map <- exp_map.df[exp_map.df$pred.cl %in% keep_cl,]
                    fg_cells <- data.frame(exp_component_name = rownames(fg_map),
                                           cl = fg_map$pred.cl)
                    
                    bg_cl <- V1.keep.cl[V1.keep.cl %in% keep_cl]
                    bg_cells <- data.frame(exp_component_name = names(bg_cl),
                                           cl = bg_cl)
                    
                    cutoff <- 0.01
                    
                    pw_de_results <- map(keep_cl, 
                                         fg_bg_pairwise_de, 
                                         fg_cells, 
                                         bg_cells)
                    
                    names(pw_de_results) <- cl.df$cluster_label[match(keep_cl, rownames(cl.df))]
                    
                    volcano_plots <- map(1:length(pw_de_results),
                                         function(x) {
                                           df <- pw_de_results[[x]]
                                           cluster_label <- names(pw_de_results)[x]
                                           de_volcano(df, 
                                                      cluster_label,
                                                      experiment,
                                                      cutoff = cutoff,
                                                      highlight_genes = immune.genes,
                                                      n_fg = sum(fg_cells$cl == keep_cl[x]),
                                                      n_bg = sum(bg_cells$cl == keep_cl[x]))
                                         })
                    
                    scatter_plots <- map(1:length(pw_de_results), 
                                         function(x) {
                                           df <- pw_de_results[[x]]
                                           cluster_label <- names(pw_de_results)[x]
                                           de_scatter(df, cluster_label, highlight_genes = immune.genes)
                                         })
                    
                    list(pw_de_results = pw_de_results,
                         volcano_plots = volcano_plots,
                         scatter_plots = scatter_plots)
                  })

names(de_results) <- experiments

map_int(de_results, function(x) length(x$volcano_plots))


walk(experiments,
     function(x) {
       
       plot_list <- de_results[[x]]$volcano_plots
        
       all_plots <- plot_grid(plotlist = plot_list,
                              nrow = 1)
       
       ggsave(paste0(x,".pdf"),
              height = 4,
              width = length(plot_list) * 4,
              useDingbats = FALSE)
       
       ggsave(paste0(x,".png"),
              height = 4,
              width = length(plot_list) * 4)
     })




#Now, let's do the same thing for the anterograde cell positive control.

#First, matching
#Let's also filter clusters based on >= 10 cells in the map.result and >= 10 cells in VISp

ant_map_table <- table(ant.map.result$map.df$pred.cl)
ant_map_cl <- ant_map_table[ant_map_table >= 10]

ant_keep_cl <- intersect(names(ant_map_cl),names(V1_keep))



#Side note: Looks like there are 17 PVM cells in the virally labeled set (pred.cl == 136)

#On to matches:

ant_fg_map <- ant.map.result$map.df[ant.map.result$map.df$pred.cl %in% ant_keep_cl,]
ant_fg_cells <- data.frame(exp_component_name = rownames(ant_fg_map),
                           cl = ant_fg_map$pred.cl)

ant_bg_cl <- V1.keep.cl[V1.keep.cl %in% ant_keep_cl]
ant_bg_cells <- data.frame(exp_component_name = names(ant_bg_cl),
                           cl = ant_bg_cl)


#Now, DE Gene computation:

ant_pw_de_results <- map(ant_keep_cl, function(cl) {
  
  cl_fg_cells <- ant_fg_cells[ant_fg_cells$cl == cl,]
  cl_bg_cells <- ant_bg_cells[ant_bg_cells$cl == cl,]
  
  cl_factor <- factor(c(paste0("fg_",cl_fg_cells$cl), 
                        paste0("bg_",cl_bg_cells$cl)))
  names(cl_factor) <- c(cl_fg_cells$exp_component_name,
                        cl_bg_cells$exp_component_name)
  
  cl_data <- cbind(ant_data[, cl_fg_cells$exp_component_name],
                   norm.dat[, cl_bg_cells$exp_component_name])
  
  cl_result <- DE_genes_pw(norm.dat = cl_data,
                           cl = cl_factor)[[1]]
  
  cl_result$gene <- rownames(cl_result)
  
  cl_result <- cl_result %>%
    arrange(padj)
  
  cl_result
  
})

names(ant_pw_de_results) <- cl.df$cluster_label[match(ant_keep_cl, rownames(cl.df))]



sum(tolower(immune.genes) %in% tolower(rownames(norm.dat)))

map(ant_pw_de_results, function(x) sum(x$padj < 0.01))


#How many overlap immune-related genes?

map(ant_pw_de_results, function(x) sum(x$padj < 0.01 & tolower(x$gene) %in% tolower(immune.genes)))


#A lot! 22 of 70 in the immune-related list.

#Plot it:

ant_volcano_plots <- map(1:length(ant_pw_de_results),
                         function(x) {
                           df <- ant_pw_de_results[[x]]
                           cluster_label <- names(ant_pw_de_results)[x]
                           
                           plot_x <- df %>%
                             mutate(color = case_when(padj < 0.01 & lfc > 0 ~ "orangered",
                                                      padj < 0.01 & lfc < 0 ~ "dodgerblue",
                                                      padj >= 0.01 ~ "black")) %>%
                             mutate(color = ifelse(padj < 0.01 & tolower(gene) %in% tolower(immune.genes),
                                                   "magenta",
                                                   color))
                           
                           plot_labels <- plot_x %>%
                             filter(padj < 0.01) %>%
                             filter(tolower(gene) %in% tolower(immune.genes))
                           
                           ggplot() +
                             geom_point(data = plot_x,
                                        aes(x = lfc,
                                            y = -log10(padj),
                                            color = color),
                                        size = 0.5) +
                             geom_text_repel(data = plot_labels,
                                             aes(x = lfc,
                                                 y = -log10(padj),
                                                 color = color,
                                                 label = gene),
                                             size = 2) +
                             scale_color_identity() +
                             xlim(-12,12) +
                             ylim(0, 100) +
                             theme_bw() +
                             ggtitle(cluster_label)
                         })

ant_scatter_plots <- map(1:length(ant_pw_de_results),
                         function(x) {
                           df <- ant_pw_de_results[[x]]
                           cluster_label <- names(ant_pw_de_results)[x]
                           
                           plot_x <- df %>%
                             mutate(color = case_when(padj < 0.01 & lfc > 0 ~ "orangered",
                                                      padj < 0.01 & lfc < 0 ~ "dodgerblue",
                                                      padj >= 0.01 ~ "black")) %>%
                             mutate(color = ifelse(padj < 0.01 & tolower(gene) %in% tolower(immune.genes),
                                                   "magenta",
                                                   color))
                           
                           plot_labels <- plot_x %>%
                             filter(padj < 0.01) %>%
                             filter(tolower(gene) %in% tolower(immune.genes))
                           
                           ggplot() +
                             geom_point(data = plot_x,
                                        aes(x = meanA,
                                            y = meanB,
                                            color = color),
                                        size = 0.5) +
                             geom_text_repel(data = plot_labels,
                                             aes(x = meanA,
                                                 y = meanB,
                                                 color = color,
                                                 label = gene),
                                             size = 2) +
                             scale_color_identity() +
                             xlim(0,12.5) +
                             ylim(0,12.5) +
                             theme_bw() +
                             ggtitle(cluster_label)
                         })


#Save the plots:

ant_plot_list <- list(ant_volcano_plots[[1]], ant_scatter_plots[[1]])

ant_all_plots <- plot_grid(plotlist = ant_plot_list,
                           ncol = 2,
                           labels = letters[11:12])

ggsave("anterograde_DE_Gene_plots_all_viruses.pdf",
       ant_all_plots,
       width = 8,
       height = 4,
       useDingbats = FALSE)

ggsave("anterograde_DE_Gene_plots_all_viruses.png",
       ant_all_plots,
       width = 8,
       height = 4)


