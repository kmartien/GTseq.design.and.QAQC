library(dplyr)
library(tidyverse)
library(swfscMisc)
source("R/functions/Compare.replicates.R")

projects <- list("RunMS51", "RunMS58")
tgt.list <- lapply(projects, function(p){
  load(paste0("results-R/", p, ".10readsMin.geno.eval.rda"))
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
replicates <- unique(LABIDs[duplicated(LABIDs)])

samp.check <-do.call('rbind',lapply(replicates, function(r){
  print(r)
  rep.tgt <- tgt[grep(substr(r, start = 3, stop = 8), tgt$Indiv),]
  locs <- unique(rep.tgt$locus)
  samp.compare <- do.call(rbind, lapply(locs, function(l){
    gts <- filter(rep.tgt, locus == l) %>% select(gt) %>% na.omit()
    return(c(locus = l, tot.genos = nrow(gts), unique.genos = nrow(unique(gts))))
  })) %>% data.frame()
  overlapping.genos <- length(which(samp.compare$tot.genos == 2))
  mismatches <- length(which(samp.compare$unique.genos > 1))
  return(c(LABID = r, overlapping.genos = overlapping.genos, mismatches = mismatches))
  #  filter(rep.tgt, locus %in% samp.compare$locus[which(samp.compare$unique.genos > 1)]) %>% arrange(locus)
})) %>% data.frame()

write.csv(samp.check, "results-raw/RunMS51.vs.RunMS58.comparisons.csv")

replicates <- filter(samp.check, mismatches > 0) %>% select(LABID)

mismatches <- do.call('rbind',lapply(replicates$LABID, function(r){
  print(r)
  rep.tgt <- tgt[grep(substr(r, start = 3, stop = 8), tgt$Indiv),]
  compare.replicates(rep.tgt)
}))

write.csv(mismatches, "results-raw/RunMS51.vs.RunMS58.mismatches.csv")
