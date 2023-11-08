library(dplyr)
library(tidyverse)
library(swfscMisc)
source("R/functions/Compare.replicates.R")

projects <- list("all.Miseq.runs", "GTseq.val", "RunMS43")
tgt.list <- lapply(projects, function(p){
  load(paste0("results-R/", p, ".10readsMin.tgt.rda"))
  tgt$Indiv <- paste0(tgt$Indiv, ".", p)
  return(tgt)
})
names(tgt.list) <- projects

tgt <- do.call(rbind, tgt.list)

unique.files <- unique(tgt$Indiv)
#MS45.ids <- unique.files[1:36]
#MS45.ids[grep("19", MS45.ids)] <- paste0("z0",MS45.ids[grep("19", MS45.ids)])
#MS45.ids[grep("18", MS45.ids)] <- paste0("z0",MS45.ids[grep("18", MS45.ids)])
#MS45.ids[grep("17", MS45.ids)] <- paste0("z0",MS45.ids[grep("17", MS45.ids)])
#MS45.ids[-grep("z", MS45.ids)] <- paste0("z00",MS45.ids[-grep("z", MS45.ids)])
#unique.files[1:36] <- MS45.ids

LABIDs <- unique.files %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]

to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- tgt[grep(substr(r, start = 3, stop = 8), tgt$Indiv),]
  mismatches <- compare.replicates(rep.tgt)
}))

write.csv(to.check, "results-raw/all.Miseq.v.GTseq.val.mismatches.csv")
