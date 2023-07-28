# This function parses the SAM file produced by BWA-MEM2 when mapping primer pairs
# to a reference genome and extracts the number of locations where each primer maps.
# The only argument is the path to a folder containing the sam files you want to
# parse. It assumes each primer pair has its own SAM file, and that they were mapped
# as paired reads.

extract.number.of.hits <- function(sam.folder.path){
  
  sam.files <- list.files(path = sam.folder.path, pattern = ".sam")
  
  primer.hits <- do.call(rbind, lapply(sam.files, function(f){
    temp <- system2(command = "grep",
                    args = c("X0", paste0(sam.folder.path, "/", f), "|", "sed"),
                    stdout = TRUE, stderr = TRUE)
    exact.matches <- do.call(cbind, lapply(1:length(temp), function(i){
      vec <- strsplit(temp[[i]], split = "\t")[[1]]
      as.numeric(strsplit(vec[grep(pattern = "X0", vec)], split = ":")[[1]][3])
    }))
    res <- data.frame(exact.matches)
    names(res) <- c("fwd", "rev")
    return(res)
  }))
  primer.names <- sapply(sam.files, function(f){
    strsplit(f, split = "[.]")[[1]][1]
  })
  primer.hits <- cbind(Locus = primer.names, primer.hits)
  names(primer.hits)[1] <- "Locus"
  
  return(primer.hits)
}