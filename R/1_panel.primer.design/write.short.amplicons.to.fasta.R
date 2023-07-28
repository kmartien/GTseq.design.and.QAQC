library(Mnov.GTseq.data)
library(dplyr)
library(strataG)

data("Mnov.gtseq.primers")
Mnov.gtseq.primers <- filter(Mnov.gtseq.primers, Orientation == "FORWARD")

ref.seqs <- read.fasta("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/data-raw/fasta/Mnov.target.amplicons.all.fasta")

# This writes all of the amplicons to a single fasta
for (i in 1:length(ref.seqs)){
  locus <- names(ref.seqs)[i]
  start.pos <- Mnov.gtseq.primers$Start[which(Mnov.gtseq.primers$Primer.name.analysis == locus)]
  stop.pos <- start.pos - 1 + Mnov.gtseq.primers$Prod.Size[which(Mnov.gtseq.primers$Primer.name.analysis == locus)]
  ref.seqs[[i]] <- ref.seqs[[i]][start.pos:stop.pos]
}
write.fasta(ref.seqs, file = "data-raw/fasta/target.amplicons.90.to.143bp.fasta")

# This writes each amplicon to its own fasta
for (i in 1:length(ref.seqs)){
  locus <- names(ref.seqs)[i]
  start.pos <- Mnov.gtseq.primers$Start[which(Mnov.gtseq.primers$Primer.name.analysis == locus)]
  stop.pos <- start.pos - 1 + Mnov.gtseq.primers$Prod.Size[which(Mnov.gtseq.primers$Primer.name.analysis == locus)]
  seq <- ref.seqs[[i]][start.pos:stop.pos]
  write.fasta(seq, file = paste0("data-raw/fasta/", locus, "-amplicon.fasta"))
}
