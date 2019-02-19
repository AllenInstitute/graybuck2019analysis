add_score_track <- function(pile_plot,
                            score_gr,
                            ucsc_loc,
                            track_ypos,
                            highlight_loc = NULL,
                            padding = c(1e5,1e5),
                            track_color = "#000000",
                            norm = "max",
                            max_val = NULL,
                            window_size = NULL,
                            window_mode = c("max","mean","median"),
                            target_color = "#B7B7B7",
                            highlight_color = "#F9ED32") {
  
  gr_target <- ucsc_loc_to_GRanges(ucsc_loc)
  target_start <- start(gr_target)
  target_end <- end(gr_target)
  
  start(gr_target) <- start(gr_target) - padding[1]
  end(gr_target) <- end(gr_target) + padding[2]
  
  # Get overlapping scores
  chr_ranges <- score_gr[as.character(seqnames(score_gr)) == as.character(seqnames(gr_target))]
  ol_ranges <- subsetByOverlaps(chr_ranges, gr_target)
  ol_len <- length(ol_ranges)
  
  if(ol_len == 0) {
    warning("No overlaps found for the score track in the target reigons")
    return(p)
  }
  
  ol_start <- start(ol_ranges)
  ol_end <- end(ol_ranges)
  ol_val <- score(ol_ranges)
  
  score_df <- map_dfr(1:length(ol_ranges), 
                      function(x) {
                        data.frame(pos = ol_start[x]:ol_end[x],
                                   val = ol_val[x])
                      })
  
  pile <- data.frame(pos = start(gr_target):end(gr_target)) %>%
    left_join(score_df)
  
  pile$val[is.na(pile$val)] <- 0
  
  if(!is.null(window_size)) {
    pile <- pile %>%
      mutate(bin = floor(pos/window_size))
    if(window_mode == "max") {
      pile <- pile %>%
        group_by(bin) %>%
        summarise(pos = bin[1]*window_size,
                  val = max(val))
    } else if(window_mode == "mean") {
      pile <- pile %>%
        group_by(bin) %>%
        summarise(pos = bin[1]*window_size,
                  val = mean(val))
    } else if(window_mode == "max") {
      pile <- pile %>%
        group_by(bin) %>%
        summarise(pos = bin[1]*window_size,
                  val = median(val))
    }
    
    pile <- pile %>%
      select(pos, val)
  }
  
  if(norm == "PM") {
    m <- sum(unlist(lapply(group_gr,length)))
    pile$val <- pile$val/m*1e6
  } else if(norm == "max") {
    pile$val <- pile$val/max(pile$val)
  }
  
  
  target_rect <- data.frame(xmin = target_start,
                            xmax = target_end,
                            ymin = track_ypos,
                            ymax = track_ypos + 1,
                            fill = target_color)
  
  chr <- as.character(seqnames(gr_target))
  
  if(!is.null(highlight_loc)) {
    hi_target <- ucsc_loc_to_GRanges(highlight_loc)
    hi_start <- start(hi_target)
    hi_end <- end(hi_target)
    
    hi_rect <- data.frame(xmin = hi_start,
                          xmax = hi_end,
                          ymin = track_ypos,
                          ymax = track_ypos + 1,
                          fill = highlight_color)
    
    pile_plot <- pile_plot +
      geom_rect(data = hi_rect,
                aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax,
                    fill = fill))
    
  }
  
  if(is.null(max_val)) {
    max_val <- max(pile$val)
  }
  
  baseline <- data.frame(x = start(gr_target),
                         xend = end(gr_target),
                         y = track_ypos,
                         yend = track_ypos,
                         color = track_color)
  
  pile_plot <- pile_plot +
    geom_segment(data = baseline,
                 aes(x = x, xend = xend,
                     y = y, yend = y,
                     color = color),
                 size = 0.1)
  
  pile$val <- track_ypos + pile$val / max_val
  pile$min <- track_ypos
  
  pile$val[pile$val > track_ypos + 1] <- track_ypos + 1
  pile_plot <- pile_plot +
    geom_ribbon(data = pile,
                aes(x = pos, ymin = min, ymax = val),
                color = NA,
                fill = track_color)
  
  
  pile_plot
  
  
}