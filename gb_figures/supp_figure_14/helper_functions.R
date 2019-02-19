fg_bg_pairwise_de <- function(cl,
                              fg_cells,
                              bg_cells) {
  
  cl_fg_cells <- fg_cells[fg_cells$cl == cl,]
  cl_bg_cells <- bg_cells[bg_cells$cl == cl,]
  
  cl_factor <- factor(c(paste0("fg_",cl_fg_cells$cl), 
                        paste0("bg_",cl_bg_cells$cl)))
  names(cl_factor) <- c(cl_fg_cells$exp_component_name,
                        cl_bg_cells$exp_component_name)
  
  cl_data <- cbind(ro_data[, cl_fg_cells$exp_component_name],
                   norm.dat[, cl_bg_cells$exp_component_name])
  
  cl_result <- DE_genes_pw(norm.dat = cl_data,
                           cl = cl_factor)[[1]]
  
  cl_result$gene <- rownames(cl_result)
  
  cl_result <- cl_result %>%
    arrange(padj)
  
  cl_result
  
}

de_volcano <- function(df, 
                       cluster_label, 
                       experiment,
                       highlight_genes,
                       cutoff = 0.01,
                       n_fg,
                       n_bg) {
  plot_x <- df %>%
    mutate(color = case_when(padj < cutoff & lfc > 0 ~ "orangered",
                             padj < cutoff & lfc < 0 ~ "dodgerblue",
                             padj >= cutoff ~ "black")) %>%
    mutate(color = ifelse(padj < cutoff & tolower(gene) %in% tolower(highlight_genes),
                          "magenta",
                          color))
  
  plot_labels <- plot_x %>%
    filter(padj < cutoff) %>%
    filter(tolower(gene) %in% tolower(highlight_genes))
  
  n_labels <- data.frame(x = c(-10, -5, -5,  5, 5, 10),
                         y = c( 90, 90, 70, 70, 90, 90),
                         label = c(n_fg, 
                                   sum(plot_x$color == "dodgerblue"),
                                   sum(plot_x$color == "magenta" & plot_x$lfc < 0),
                                   sum(plot_x$color == "magenta" & plot_x$lfc > 0),
                                   sum(plot_x$color == "orangered"),
                                   n_bg),
                         color = c("black","dodgerblue","magenta","magenta","orangered","black"))
  
  ggplot() +
    geom_text(data = n_labels,
              aes(x = x, y = y,
                  color = color,
                  label = label)) +
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
                    size = 4) +
    scale_color_identity() +
    scale_x_continuous("log2(Fold Change)", limits = c(-12, 12)) +
    scale_y_continuous("log10(p-value)", limits = c(0,100)) +
    theme_bw(8) +
    ggtitle(paste(experiment,cluster_label))
}

de_scatter <- function(df, 
                       cluster_label,
                       highlight_genes) {

  plot_x <- df %>%
    mutate(color = case_when(padj < 0.01 & lfc > 0 ~ "orangered",
                             padj < 0.01 & lfc < 0 ~ "dodgerblue",
                             padj >= 0.01 ~ "black")) %>%
    mutate(color = ifelse(padj < 0.01 & tolower(gene) %in% tolower(highlight_genes),
                          "magenta",
                          color))
  
  plot_labels <- plot_x %>%
    filter(padj < 0.01) %>%
    filter(tolower(gene) %in% tolower(highlight_genes))
  
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
}