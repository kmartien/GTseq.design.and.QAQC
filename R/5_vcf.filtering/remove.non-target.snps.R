library(vcfR)
library(vcftoolsR)
library(tidyverse)
library(dplyr)
library(Mnov.GTseq.data)
source("R/vcf2genos/0 Load vcf2genos functions.R")

# description should be the name of the vcf file (without extension) that you want to filter
description = "ddRAD.528locs"

# If you don't have a bed file in a data package, you can read one in from a .csv
data("SWFSC.bed")
bed <- SWFSC.bed
bed <- read.csv("data-raw/528locs.allSNPs.bed.csv")

vcf <- read.vcfR(paste0("vcf/", description, ".vcf"))

tgt <- create.tgt(vcf)
loci <- unique(tgt$locus)

sites.to.keep <- paste0(bed$locus, "_", bed$stop)
non.target.loci <- filter(tgt, !locus %in% sites.to.keep) %>% select(c(locus, CHROM, POS)) %>%
  distinct()

write.table(non.target.loci[,c(2,3)], file="vcf/temp/non.target.locs.to.remove2.txt", 
            row.names = FALSE, col.names=FALSE, quote = FALSE)
vcftools.rmPOS(paste0("vcf/", description), 
                 paste0("vcf/", description, ".allSNPs"),
                 "vcf/temp/non.target.locs.to.remove2.txt")
#Use this line if you don't have vcftools installed in a way that can be called from R
#paste0("vcftools --vcf vcf/", description, ".vcf --out vcf/", description, ".targetSNPs --exclude-positions vcf/temp/non.target.locs.to.remove2.txt --recode --recode-INFO-all")

fname <- paste0(description, ".allSNPs.recode")
vcf <- read.vcfR(paste0("vcf/", fname, ".vcf"))
