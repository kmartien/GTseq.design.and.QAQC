library(tidyverse)
library(dplyr)
source("R/functions/mplot2tgt.R")
source("R/functions/Compare.replicates.R")

project <- "GTseq.val"

AB.min.het <- 3/7
AB.max.homo <- 2/8
min.AR.het <- 3/10
max.AR.homo <- 2/10
min.read.depth <- 10
num.locs <- 368
min.genos.per.ind <- num.locs * 0.8
why.removed <- list()

tgt <- mplot2tgt(project = project, AB.min.het = AB.min.het, AB.max.homo = AB.max.homo,
                 min.read.depth = min.read.depth)

# compare replicates
LABIDs <- unique(tgt$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]
mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- filter(tgt, Indiv %in% c(r,paste0(r,"b")))
  mismatches <- compare.replicates(rep.tgt)
}))
if(nrow(to.check > 0)) {
  print("Some replicates have mismatched genotypes")
  print(paste0("Mismatches saved to results-R/", project, ".genotype.mismatches.rda"))
  save(mismatches.to.check, file = paste0("results-R/", project, ".genotype.mismatches.rda"))
}

# Identify samples with Ns or Xs in their haplotypes or more than 2 haplotypes
questionable.hap <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("N", tgt$gt[i],)) > 0 || length(grep("X", tgt$gt[i],)) > 0 
         || tgt$num.haps[i] > 2, TRUE, FALSE) 
})
genos.to.check <- filter(tgt, questionable.hap == TRUE)
if(nrow(to.check > 0)) {
  print("Some samples with Ns or Xs in their haplotypes or more than 2 haplotypes")
  print(paste0("Questionable genotypes saved to results-R/", project, ".genos.to.check.rda"))
  save(genos.to.check, file = paste0("results-R/", project, ".genos.to.check.rda"))
}

# summarize individual data
missing.data.ind <- data.frame(table(tgt$Indiv[!is.na(tgt$gt)])) %>%
  mutate(missing = num.locs-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$genos >= min.genos.per.ind))
write.csv(missing.data.ind, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.num.genos.per.ind.csv"))
num.inds <- nrow(missing.data.ind)

# summarize locus data
missing.data.loc <- data.frame(table(tgt$locus[!is.na(tgt$gt)])) %>%
  mutate(missing = num.inds-Freq)
names(missing.data.loc)[1] <- "locus"
length(which(missing.data.loc$missing > round(num.inds * 0.2)))

geno.table <- tgt.2.geno.table(tgt)

geno.table$num.genos <- do.call(rbind, lapply(1:nrow(geno.table), function(i){
  length(which(!is.na(geno.table[i,2:ncol(geno.table)])))
}))

genos.per.locus <- do.call(rbind, lapply(2:(ncol(geno.table)-1), function(l){
  length(which(!is.na(geno.table[1:nrow(geno.table),l])))
}))
genos.per.locus <- distinct(data.frame(colnames(geno.table)[2:(ncol(geno.table)-1)], genos.per.locus))
names(genos.per.locus) <- c("locus", "genos")

loci <- unique(genos.per.locus$locus)
genotypes.per.locus <- data.frame(loci, do.call(rbind, lapply(loci, function(loc){
  haps <- filter(tgt, locus == loc) %>% filter(!is.na(gt))
  length(unique(haps$gt))
})))
names(genotypes.per.locus) <- c("locus", "num.haps")
loc.sum <- left_join(genotypes.per.locus, genos.per.locus)
write.csv(loc.sum, file = paste0("results-raw/", project, ".locus.summary.csv"))
save(loc.sum, file = paste0("results-R/", project, ".locus.summary.rda"))

save(geno.table, tgt, file = paste0("results-R/", project, ".", min.read.depth, "readsMin.tgt.rda"))

