library(dplyr)
library(tidyverse)
library(Mnov.GTseq.data)

run <- "all"
project <- paste0(run,".10readsMin")
load(paste0("results-R/", project, ".tgt.rda"))
data("SWFSC.bed")
num.locs <- length(unique(tgt$locus))

locs.not.in.tgt <- data.frame(SWFSC.bed$locus[-which(SWFSC.bed$locus %in% tgt$locus)])

# calculate missingness by individual
missing.data.ind <- data.frame(table(tgt$Indiv)) %>%
  mutate(missing = num.locs-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")

# remove individuals with 50 or fewer genotypes
rejected.inds <- missing.data.ind$labID[which(missing.data.ind$genos <= 50)]
tgt <- tgt[-which(tgt$Indiv %in% rejected.inds),]
num.inds <- length(unique(tgt$Indiv))

# calculate missingness by locus
missing.data.loc <- data.frame(table(tgt$locus)) %>%
  mutate(missing = num.inds-Freq)
names(missing.data.loc)[1] <- "locus"
length(which(missing.data.loc$Freq < (num.inds * 0.75)))

# remove loci genotyped in less than 75% of the individuals
rejected.locs <- missing.data.loc$locus[which(missing.data.loc$Freq < (num.inds * 0.75))]
tgt <- tgt[-which(tgt$locus %in% rejected.locs),]
num.locs <- length(unique(tgt$locus))

