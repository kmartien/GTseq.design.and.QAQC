library(dplyr)
library(vcfR)
library(tidyverse)
library(pegas)
source("R/functions/define_snp_windows.R")
source("R/functions/vcf2allelefreqs.R")

vcf <- read.vcfR("vcf/Mnov_filtered_SNPs.vcf", convertNA = T)

## identify windows
window.list <- define_snp_windows(tidy.vcf, allele.freqs)

## collapse list into a single dataframe
all.windows <- do.call(rbind, window.list)


