# This function determines the number of reads in fastq.gz files. It counts the 
# number of lines in each file and divides by four (since there are four lines 
# per read in a fastq). It can be used with either single or paired reads. For 
# paired reads it counts the number of reads in the forward file.

count.reads.in.fastq <- function(fastq.dir, fastq.files){ #, read.type){
  
    tot.reads <- do.call(rbind, lapply(1:length(fastq.files), function(f){
      cnt <- system2(command = "zcat",
                     args = c("<", paste0(fastq.dir, "/", fastq.files[f]), "|", "wc", "-l"),
                     stdout = TRUE, stderr = TRUE)
      return(as.numeric(cnt)/4)
    }))
    
  return(data.frame(Sample = fastq.files, tot.reads = tot.reads))
}