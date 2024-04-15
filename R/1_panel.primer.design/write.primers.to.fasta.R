library(dplyr)
#library(Mnov.GTseq.data)

#data("Mnov.gtseq.primers")
#primers <- Mnov.gtseq.primers
fwd.primers <- filter(selected_rows, Orientation == "FORWARD") %>% arrange(seq_id)
rev.primers <- filter(selected_rows, Orientation == "REVERSE") %>% arrange(seq_id)

for (i in 1:nrow(fwd.primers)){
  cat(paste0("> ", fwd.primers$locuschrom[i],".", fwd.primers$ampstart[i],"-", fwd.primers$ampend[i],"\n"), file = paste0("JA_DATA/fasta/primers/",fwd.primers$locuschrom[i],".", fwd.primers$ampstart[i],"-", fwd.primers$ampend[i],".R1.fasta"))
  cat(paste0(fwd.primers$Seq[i], "\n"), file = paste0("JA_DATA/fasta/primers/",fwd.primers$locuschrom[i],".", fwd.primers$ampstart[i],"-", fwd.primers$ampend[i],".R1.fasta"),append = TRUE)
}

for (i in 1:nrow(fwd.primers)){
  cat(paste0("> ", fwd.primers$locuschrom[i],".", fwd.primers$ampstart[i],"-", fwd.primers$ampend[i],"\n"), file = paste0("JA_DATA/fasta/primers/",fwd.primers$locuschrom[i],".", fwd.primers$ampstart[i],"-", fwd.primers$ampend[i],".R2.fasta"))
  cat(paste0(rev.primers$Seq[i], "\n"), file = paste0("JA_DATA/fasta/primers/",fwd.primers$locuschrom[i],".", fwd.primers$ampstart[i],"-", fwd.primers$ampend[i],".R2.fasta"),append = TRUE)
}
