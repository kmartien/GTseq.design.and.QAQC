library(dplyr)
library(Mnov.GTseq.data)

data("Mnov.gtseq.primers")
primers <- Mnov.gtseq.primers
fwd.primers <- filter(primers, Orientation == "FORWARD") %>% arrange(Primer.name.analysis)
rev.primers <- filter(primers, Orientation == "REVERSE") %>% arrange(Primer.name.analysis)

for (i in 1:nrow(fwd.primers)){
  cat(paste0("> ", fwd.primers$Primer.name.analysis[i], ".", fwd.primers$locus[i], "\n"), file = paste0("data-raw/fasta/primers/",fwd.primers$Primer.name.analysis[i],".R1.fasta"))
  cat(paste0(fwd.primers$loc.Seq[i], "\n"), file = paste0("data-raw/fasta/primers/",fwd.primers$Primer.name.analysis[i],".R1.fasta"), append = TRUE)
}

for (i in 1:nrow(fwd.primers)){
  cat(paste0("> ", rev.primers$Primer.name.analysis[i], ".", rev.primers$locus[i], "\n"), file = paste0("data-raw/fasta/primers/",fwd.primers$Primer.name.analysis[i],".R2.fasta"))
  cat(paste0(rev.primers$loc.Seq[i], "\n"), file = paste0("data-raw/fasta/primers/",fwd.primers$Primer.name.analysis[i],".R2.fasta"), append = TRUE)
}
