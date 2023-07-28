library(dplyr)
library(Mnov.GTseq.data)

#data("AS.gtseq")
#data("ind.info")
#data("Mnov.strata")
#data("gtseq.prod.summary")
#data("gtseq.val.genos")
#data("tissue.archive.data")
data("GTseq.samps.final")
#data("id.key")

load("results/GTseq.sample.success.summary.rda")

successes <- filter(read.sum, success == TRUE)
inds <- left_join(successes, GTseq.samps.final)
table(inds$GTseq.stratum)
