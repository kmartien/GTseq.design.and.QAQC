remove.rejected.loci <- function(tgt, locus.ARs){
  rejected.locs <- locus.ARs$locus[which(locus.ARs$status == "Reject")]
  temp <- tgt[-which(tgt$locus %in% rejected.locs),]
  return(temp)
}