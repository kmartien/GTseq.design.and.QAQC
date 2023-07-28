library(dplyr)
library(Mnov.GTseq.data)
load("R/functions/extract.number.of.hits.R")

data("Mnov.GTSEEK.panel")

sam.folder.path <- "data-raw/sam.files/Paired.primers.mapped.to.genome"

primer.matches <- extract.number.of.hits(sam.folder.path)

primer.matches$in.panel <- primer.matches$Locus %in% Mnov.GTSEEK.panel$short

GTSEEK.primer.test <- read.csv("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/SNPs/GTSEEK_validation_run/HB-Whale_PrimerTest.csv")
primer.matches <- left_join(primer.matches, GTSEEK.primer.test)

save(primer.matches, file = "data/primer.num.matches.to.genome.rda")
write.csv(primer.matches, file = "results-raw/primer.num.matches.to.genome.csv")
