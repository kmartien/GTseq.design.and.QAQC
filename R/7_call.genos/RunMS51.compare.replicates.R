library(dplyr)
library(tidyverse)

project <- "RunMS51.allSamps.10readsMin"
load(paste0("results-R/", project, ".tgt.rda"))

LABIDs <- unique(tgt$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]

to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- filter(tgt, Indiv %in% c(r,paste0(r,"b")))
  mismatches <- compare.replicates(rep.tgt)
}))
