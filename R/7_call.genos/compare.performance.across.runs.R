library(tidyverse)
library(dplyr)
library(Mnov.GTseq.data)

min.reads <- 10

run.names <- c("GTseq.val", "GTseq.prod", "RunMS43", "RunMS45.90s", "RunMS51", "RunMS54", "RunMS58")

all.runs <- lapply(run.names, function(r){
  load(paste0("results-R/", r, ".", min.reads, "readsMin.geno.eval.rda"))
  genos.per.ind <- data.frame(table(tgt$Indiv[!is.na(tgt$gt)]))
  names(genos.per.ind) <- c("LABID", "genos")
  return(genos.per.ind)
})

all.samples <- do.call(rbind, lapply(all.runs, function(r){
  return(r)
})) %>% select(LABID) %>% distinct(LABID)

for (r in 1:length(all.runs)){
  all.samples <- left_join(all.samples, all.runs[[r]], by = "LABID")
}
names(all.samples)[2:8] <- run.names
write.csv(all.samples, file = "results-raw/sample.performance.across.all.runs.csv")
