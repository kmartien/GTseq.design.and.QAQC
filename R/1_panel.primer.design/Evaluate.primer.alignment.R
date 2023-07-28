library(dplyr)
library(Mnov.GTseq.data)

data("Mnov.GTSEEK.panel")

sam.folder <- "Paired.primers.mapped.to.genome"

sam.files <- list.files(path = paste0("data-raw/sam.files/", sam.folder), pattern = ".sam")

primer.matches <- do.call(rbind, lapply(sam.files, function(f){
  temp <- system2(command = "grep",
                 args = c("X0", paste0("data-raw/sam.files/", sam.folder, "/", f), "|", "sed"),
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
primer.matches <- cbind(Locus = primer.names, primer.matches)
names(primer.matches)[1] <- "Locus"
primer.matches$in.panel <- primer.matches$Locus %in% Mnov.GTSEEK.panel$short

GTSEEK.primer.test <- read.csv("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/SNPs/GTSEEK_validation_run/HB-Whale_PrimerTest.csv")
primer.matches <- left_join(primer.matches, GTSEEK.primer.test)

save(primer.matches, file = "data/primer.num.matches.to.genome.rda")
write.csv(primer.matches, file = "results-raw/primer.num.matches.to.genome.csv")
