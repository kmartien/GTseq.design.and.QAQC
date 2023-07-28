geno.smry.by.locus <- function(tgt, fn.prefix){
  
  tgt <- tgt %>%
    select(CHROM, locus, Indiv, REF, ALT, gt) %>%
    separate(gt, c("haplo.x", "haplo.y"), sep = "/", remove = FALSE)
  
  tgt$type <- sapply(1:nrow(tgt), function(i){
    ifelse(tgt$haplo.x[i] != tgt$haplo.y[i], "HET", 
           ifelse(tgt$haplo.x[i] == tgt$REF[i], "R.HOMO", 
                  ifelse(tgt$haplo.x[i] == tgt$ALT[i], "A.HOMO", "NA")))
  })
  
  loci <- unique(tgt$locus)
  
  possible.types <- c("R.HOMO", "HET", "A.HOMO", "NA")
  
  loc.sum <- data.frame(do.call(rbind, lapply(loci, function(i){
    tgt.i <- filter(tgt, locus == i)
    sapply(possible.types, function(t){
      length(which(tgt.i$type == t))
    })
  })))
  
  loc.sum <- cbind(loci, loc.sum)
  names(loc.sum)[1] <- "locus"
  
  write.csv(loc.sum, file = paste0(fn.prefix, "_geno.smry.by.locus.csv"))
  return(loc.sum)
}
