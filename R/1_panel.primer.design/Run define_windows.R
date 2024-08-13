library(dplyr)
library(vcfR)
library(tidyverse)
library(pegas)
source("R/functions/define_snp_windows.R")
source("R/functions/vcf2allelefreqs.R")
source("R/functions/evaluate_window_clusters.R")
source("R/functions/read.chosen.loc.seqs.from.reference.R")
source("R/functions/extract.target.fragments.batchprimer3.R")

vcf <- read.vcfR("vcf/Mnov_filtered_SNPs.vcf", convertNA = T)

# convert vcf to a tibble and add a column 'locus' that combines CHROM and POS
# of the SNPs
tidy.vcf <- vcfR2tidy(vcf, single_frame = TRUE, 
                      info_fields = c("DP"), format_fields = c("GT"))$dat %>%
  mutate(locus = paste(CHROM, POS, sep= "_")) %>%
  relocate(locus, .after = POS)

## use these lines to limit length for testing purposes
load("define.windows.locs.2.exclude.rda")
tidy.vcf <- tidy.vcf[1:100000,] %>% filter(!locus %in% locs.2.exclude$locus)

# I need the following line to get rid of the '.1' that appears at the end of my
# CHROM names in the vcf. You may not need it. CHROM in tidy.vcf must match
# CHROM in reference genome.
tidy.vcf$CHROM <- substr(tidy.vcf$CHROM, 1, 12)

## identify windows with define_snp_windows(vcf, window.sz = 134, min.SNPs = 1, max.SNPs = 5)
# This function identifies all possible 134bp windows containing between 1 and 5
# SNPs (those are default values) by checking for SNPs within 134 downstream of 
# each focal SNP. Windows are identified by their start and stop positions. 
# Note that windows can overlap. For instance, if there are SNPs at positions 13, 
# 56, 104, and 160, there will be 4 windows: 13-104, 56-160, 104-160, and 160-160.
window.list <- define_snp_windows(tidy.vcf)

## evaluate overlapping windows with evaluate_window_clusters(window.list)
# This script identifies 'clusters' as being groups of windows where the stop
# position of one windows is after the start position of the next one. It then 
# one window from each cluster to keep based on the mean heterozygosit and number
# of SNPs in each window and returns one window per cluster.
final.windows <- do.call(rbind, evaluate_window_clusters(window.list))

chosen.loc.seqs.file <- "data/chosen.loc.seqs.rda" # Where do you want to store chosen.loc.seqs?
ref.genome.location <- "/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1."

# read reference genome and only keep the bits that contain your final windows
if (file.exists(chosen.loc.seqs.file)) load(chosen.loc.seqs.file) else {
  chosen.locs.seqs <- read.chosen.loc.seqs.from.reference(final.windows, ref.genome.location)
}

# These next two lines are to correct for the fact that I'm only processing a subset
# of my data to make debugging go faster. You can skip them
chosen.loc.seqs2 <- chosen.loc.seqs[which(names(chosen.loc.seqs) %in% final.windows$CHROM)]
final.windows2 <- filter(final.windows, CHROM %in% names(chosen.loc.seqs2))

# Use function extract.target.fragments(tidy.vcf, chosen.loc.seqs) to extract
# sequence fragments from the chosen loci to send to Primer 3 for primer design.
# Each fragment is 259bp long and is centered on the center of a SNP window. 
# Square brackets are added before the first SNP in the window and after the last
# one to identify for Primer3 the portion of the sequence that should be 
# amplified by the primers it designs.
target.frags <- extract.target.fragments(tidy.vcf, chosen.loc.seqs)

save(tidy.vcf, vcf, window.list, final.windows, target.frags, file = "data/selected.windows.data.rda")
# If you've already run the above code and want to skip it, just load these files
load("data/selected.windows.data.rda")
load("data/chosen.loc.seqs.rda")

# Write target amplicons to files.  BatchPrimer3 only accepts 500 amplicons at a
# time, so the amplicons are written to multiple files in groups of 500
# NOTE: this code uses append, so if you run it multiple times you'll just keep
# adding stuff to the end of the existing files rather than overwriting them
lapply(1:length(target.frags), function(i){
  fname = paste("data-raw/Mnov_target_amplicons_",(floor(i/500)+1),".fasta",sep="")
  write(paste(">",names(target.frags)[i],sep=" "),file=fname,append = TRUE)
  write(paste0(target.frags[[i]], collapse = ""), file = fname, append = TRUE)
  names(target.frags)[i]
})
  
