microhaplot.2.tgt <- function(genos){

  genos <- mutate(genos, id.loc = paste0(id, "-", locus))
  genos$to.remove <- FALSE
  
  genos$num.haps <- sapply(1:nrow(genos), function(i){
    filter(genos, id.loc == genos$id.loc[i]) %>% 
      select(rank) %>% max()
  })
  
  genos$to.remove[which(genos$rank > 4)] <- TRUE
  genos.filtered <- filter(genos, to.remove == FALSE) 
  
  # create tgt-like data structure
  tgt <- data.frame(do.call(rbind, lapply(unique(genos.filtered$id.loc), function(x){
    df <- filter(genos.filtered, id.loc == x) %>% arrange(rank)
    loc <- df$locus[1]
    ind <- df$id[1]
    num.haps <- df$num.haps[1]
    haps <- df$haplo
    if(length(haps) < 4) haps <- c(haps, rep(NA, (4-length(haps))))
    depth <- df$depth
    if(length(depth) < 4) depth <- c(depth, rep(0, (4-length(depth))))
    # If a replicated sample is homozygous, microhaplot treats the calls from the differnt
    # sam.files as different haplotypes; this fixes that
#    if(haps[1] == haps[2]){
#      depth[1] <- depth[1] + depth[2]
#      depth[2] <- 0
#      haps[2] <- NA
#    }
    res <- c(loc, ind, haps, depth, num.haps)
    names(res) <- c("locus", "Indiv", "haplo.1", "haplo.2", "haplo.3", "haplo.4", "gt_RO", "A2", "A3", "A4", "num.haps.called")
    return(res)
  })))
  
  tgt$gt_RO <- as.integer(tgt$gt_RO)
  tgt$A2 <- as.integer(tgt$A2)
  tgt$A3 <- as.integer(tgt$A3)
  tgt$A4 <- as.integer(tgt$A4)
  
  tgt$ratio.x <- tgt$gt_RO/(tgt$gt_RO + tgt$A2)
  tgt$ratio.y <- tgt$A2/(tgt$gt_RO + tgt$A2)
  tgt$ratio.x <- as.numeric(round(tgt$ratio.x, 2))
  tgt$ratio.y <- as.numeric(round(tgt$ratio.y, 2))
  
  tgt$gt <-  ifelse(is.na(tgt$haplo.2), paste(tgt$haplo.1, tgt$haplo.1, sep="/"), paste(tgt$haplo.1, tgt$haplo.2, sep="/"))
  tgt$depth.x <- as.numeric(tgt$gt_RO)
  tgt$depth.y <- as.numeric(tgt$A2) 
  
  tgt$REF <- tgt$ALT <- tgt$QUAL <- NA

  return(tgt)
}
