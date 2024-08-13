evaluate_window_clusters <- function(window.list){
  
  lapply(window.list, function(ch){
    ch$cluster <- 1
    if(nrow(ch) == 1) return(ch)
    for (i in 2:nrow(ch)){
      if (ch$start[i] <= ch$stop[i-1]) {
        ch$cluster[i] <- ch$cluster[i-1]
        } else ch$cluster[i] <- ch$cluster[i-1]+1
    }
    # Choosing the window with the highest theta isn't the right choice. In some cases
    # there are multiple windows that include a high_h SNP and some lower h SNPs.
    # In these cases, the window with the fewest SNPs has the highest h since the
    # high_h SNP gets more weight, but I think it would be better to include the
    # window with the most SNPs that includes the high_h SNP. Here's an example:
    
#    window                 start   stop    num.snps mean_h     cluster
#    RYZJ01000004.1_5762746 5762746 5762791 3        0.24133504 15
#    RYZJ01000004.1_5762754 5762754 5762791 2        0.24853246 15
    
    # Here are the allele frequencies of the 3 SNPs involved:
#    Locus                   allele.1  allele.2
#    RYZJ01000004.1_5762746    36        8
#    RYZJ01000004.1_5762754    31        13
#    RYZJ01000004.1_5762791    38        6
    
    # Instead find the window with the highest mean_h, get it's stop position, filter
    # for all windows with the same stop position, and choose the one with the lowest start.
    # That gets the window with the most SNPs that includes the window with the 
    # highest mean_h
    
    chosen.windows <- do.call(rbind, lapply(1:length(unique(ch$cluster)), function(clust){
      windows.in.cluster <- filter(ch, cluster == clust)
      max.mean.theta <- max(windows.in.cluster$mean.theta_h)
      # select window(s) with highest mean theta (there may be ties)
      remaining.windows <- filter(windows.in.cluster, mean.theta_h == max.mean.theta)
      most.SNPs <- max(remaining.windows$n.snps)
      # of the windows with the higest theta, choose the one with the most SNPs, 
      # then arrange by amplicon length and take the shortest one (shorter amplicons are
      # more likely to have primers successfully designed)
      remaining.windows <- filter(remaining.windows, n.snps == most.SNPs) %>% arrange(stop-start)
      return(remaining.windows[1,])
    }))
    chosen.windows$locus <- paste(chosen.windows$CHROM, chosen.windows$start, chosen.windows$stop, sep = "_")
    return(chosen.windows)
  })
}