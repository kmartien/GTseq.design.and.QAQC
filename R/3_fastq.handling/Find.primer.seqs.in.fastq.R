# This scripts counts the number of reads in each fastq.gz file that contain
# each forward and reverse locus primer

library(spgs)

primer.seqs.in.fastq <- function(fastq.dir, fastq.files, primers, read.type) {
  
  all.matches <- do.call(rbind, lapply(1:length(fastq.files), function(f){
    fname <- fastq.files[f]
    print(paste0(f, " ", fname))
    if(read.type == "paired"){  #for paired end reads, search for the reverse primer in the reverse file
      fname.rev <- gsub("_R1", "_R2", fname)
    } else fname.rev <- fname
    matches <- do.call(rbind, lapply(1:nrow(primers), function(p){
      fwd.seq <- primers$fwd[p]
      fwd.count <- system2(command = "zgrep",
                           args = c(fwd.seq, paste0(fastq.dir, "/", fname), "-o", "|", "wc", "-l"),
                           stdout = TRUE, stderr = TRUE)
      if(read.type == "paired"){
        rev.seq <- primers$rev[p]
      } else rev.seq <- reverseComplement(primers$rev[p], case = "as is")
      rev.count <- system2(command = "zgrep",
                           args = c(rev.seq, paste0(fastq.dir, "/", fname.rev), "-o", "|", "wc", "-l"),
                           stdout = TRUE, stderr = TRUE)
      return(c(fwd = as.numeric(fwd.count), rev = as.numeric(rev.count)))
    }))
    matches <- bind_cols(locus = primers$locus, Sample = fname, data.frame(matches))
    return(matches)
  }))
  return(all.matches)
}
