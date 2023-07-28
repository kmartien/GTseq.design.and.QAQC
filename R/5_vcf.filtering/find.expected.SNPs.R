library(ape)
library(dplyr)

GTseq.loci <- read.csv("data-raw/rename_chrom.csv", header = FALSE)
names(GTseq.loci) <- c("locus", "new.name")
#GTseq.loci <- separate_wider_delim(GTseq.loci, cols = locus, delim = "_", names = c("CHROM", "POS"), cols_remove = FALSE)
#GTseq.loci$POS <- as.numeric(GTseq.loci$POS)

seqs <- read.FASTA("data-raw/fasta/Mnov.GTseq.528.allSNPsasNs.fasta")

var.sites <- lapply(seqs, function(i){
  which(i == "f0")
})

cbind("x", var.sites[[1]])

bed.file <- data.frame(do.call(rbind, lapply(1:length(var.sites), function(i){
  cbind(names(var.sites)[i], var.sites[[i]])
})))
names(bed.file) <- c("locus","stop")
bed.file$stop <- as.numeric(bed.file$stop)
bed.file$start <- bed.file$stop-1

bed.file <- left_join(bed.file, GTseq.loci) %>% select(c(new.name, start, stop))
names(bed.file)[1] <- "locus"

write.csv(bed.file, file = "data-raw/528locs.allSNPs.bed.csv", row.names = FALSE)
