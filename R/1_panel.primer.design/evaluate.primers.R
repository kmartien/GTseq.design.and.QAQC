library(tidyverse)
library(seqinr)

sam.file <- "data-raw/bam.files/Bphy.primers.mapped.to.genome.sam"
load(file = "data-raw/primer3/Bphy.final.windows.and.primers.rda")
amplicons <- read.fasta("data-raw/fasta/Bphy.amplicons.fasta", whole.header = TRUE)

sam <- readLines(sam.file)
sam <- sam[-grep("@", sam)] #remove the header

hits <- data.frame(do.call(rbind, lapply(1:length(sam), function(l){
  x <- strsplit(sam[l], split = '\t')[[1]]
  best.hits <- x[grep("X0", x)]
  best.hits <- strsplit(best.hits, split = ":")[[1]][3]
  s.h <- grep("X1", x)
  if (length(s.h) == 0) suboptimal.hits <- NA else 
    suboptimal.hits <- strsplit(x[grep("X1", x)], split = ":")[[1]][3]
  return(c(SEQUENCE =  x[10], best.hits = best.hits, suboptimal.hits = suboptimal.hits))
})))

ns.per.amplicon <- data.frame(do.call(rbind, lapply(amplicons, function(amp){
  return(c(names(amp), ns = length(grep("n", amp))))
}))) %>% rownames_to_column(var = "locus")

primer.summary <- left_join(select(final.primers, c(locus, orientation, SEQUENCE, PRODUCT_SIZE, length)),
                            hits) %>% left_join(select(filtered.windows, c(locus, mean.theta_h, n.snps))) %>% left_join(ns.per.amplicon)
