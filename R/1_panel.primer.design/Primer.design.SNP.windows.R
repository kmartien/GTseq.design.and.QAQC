library(seqinr)
library(vcfR)
library(tidyr)
library(dplyr)
library(tidyverse)
library(pegas)
source("R/functions/define_snp_windows.R")
source("R/functions/vcf2allelefreqs.R")
source("R/functions/evaluate_window_clusters.R")
source("R/functions/read.chosen.loc.seqs.from.reference.R")
source("R/functions/get.primer3.primers.R")

description <- "Bphy"

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
tidy.vcf$CHROM <- do.call('rbind',strsplit(tidy.vcf$CHROM,".",fixed=TRUE))[,1]

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

chosen.loc.seqs.file <- "data/Fernandas.chosen.loc.seqs.rda" # Where do you want to store chosen.loc.seqs?
ref.genome.location <- "/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.fna"
ref.genome.location <- "data-raw/GCA_023338255.1_SBiKF_Bphy_ph2_genomic.fna"

# read reference genome and only keep the bits that contain your final windows
ref.genome <- read.fasta(ref.genome.location)
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs <- ref.genome[which(names(ref.genome) %in% final.windows$CHROM)]
save(chosen.loc.seqs,file=chosen.loc.seqs.file)
rm(ref.genome)
load("data/Fernandas.chosen.loc.seqs.rda")

# These next two lines are to correct for the fact that I'm only processing a subset
# of my data to make debugging go faster. You can skip them
chosen.loc.seqs2 <- chosen.loc.seqs[which(names(chosen.loc.seqs) %in% final.windows$CHROM)]
final.windows2 <- filter(final.windows, CHROM %in% names(chosen.loc.seqs2))

###########################################################################
#
#         The following lines require use of the commandline version
#         of Primer3. You can install it using Homebrew. If you don't 
#         already have Homebrew installed, you can install it here:
#         
#         https://brew.sh/
#         
#         Then open a Terminal window and use this command to install Primer3:
#         
#         $ brew install primer3
#         
#         You can find this command and info about installing Primer3 at
#         this site: https://formulae.brew.sh/formula/primer3
#
#         When I installed Primer3, the folder containing the configuration
#         files (primer3_config/) didn't end up in a place where Primer3
#         could find it. I had to download the source code for Primer3, find
#         the primer3_config folder (for me it was in primer3-2.4.0/src/), and
#         put it in a location of my choice (config.folder). I then pass 
#         that location to Primer3 in my script. It's clunky, but it works.
#         
###########################################################################


# Use function extract.target.fragments.batchprimer3 (tidy.vcf, final.windows, 
# chosen.loc.seqs) to extract sequence fragments from the chosen loci, creates 
# an input file for Primer3, and calls Primer3 for primer design. Each fragment
# is 259bp long and is centered on the center of a SNP window. The portion of 
# the sequence that should be amplified by the primers Primer3 designs is 
# specified for each window in the Primer3 input file. The function requires as 
# an argument the location of the Primer3 configuration folder. The function
# captures the output from Primer3 and reformats it into the dataframe primer.deets.
# It also saves the sequences of the target fragments that were used in primer
# design, which are needed further down the script in order to extract the 
# amplicon of each primer pair
config.folder = "/Users/Shared/KKMDocuments/Documents/Github.Repos/primer3_config/"
primer3.input.fname = 
  paste("data-raw/primer3/", description, "_primer3.input_",Sys.Date(),".txt",sep="")

primer.res <- get.primer3.primers(description, tidy.vcf, final.windows, 
                            chosen.loc.seqs, config.folder, primer3.input.fname)
save(primer.res, file = paste0("data-raw/primer3/", description, ".primer.results.", Sys.Date(), ".rda"))
primer.deets <- primer.res$primer.deets 
targets <- data.frame(primer.res$targets)
names(targets) <- c("locus", "fragment")

#count the number of primer pairs successfully designed for each locus.  Remove loci that don't have primers
final.windows <- group_by(primer.deets, locus) %>% 
  summarise(n.primers = length(orientation)/2) %>%
  right_join(final.windows) %>% filter(!is.na(n.primers))

# Calculate distance between adjacent windows. Windows that don't have another 
# window within 100,000 bp are good. For those that do, decide which windows to keep
final.windows.bkup <- final.windows
final.windows$status <- NA

for(i in 1:nrow(final.windows)){
  close.snps <- filter(final.windows, CHROM == final.windows$CHROM[i]) %>%
    filter(start > (final.windows$start[i] - 100000)) %>%
    filter(start < (final.windows$start[i] + 100000))
  if(nrow(close.snps) == 1) final.windows$status[i] <- "keep" else {
    to.keep <- filter(close.snps, n.snps == max(close.snps$n.snps))
    to.keep <- filter(to.keep, mean.theta_h == max(to.keep$mean.theta_h)) %>%
      select(locus)
    to.ditch <- close.snps$locus[which(close.snps$locus != to.keep$locus)]
    final.windows$status[which(final.windows$locus == to.keep$locus)] <- "keep"
    final.windows$status[which(final.windows$locus %in% to.ditch)] <- "ditch"
  }
}

filtered.windows <- filter(final.windows, status == "keep")
final.primers <- filter(primer.deets, locus %in% filtered.windows$locus) %>%
  filter(index == 0)
save(filtered.windows, final.primers, file = "data-raw/primer3/Bphy.final.windows.and.primers.rda")

# write primer sequences to fasta file
primer.fasta.fname <- paste0("data-raw/fasta/", description, ".primers.fasta")
fwd.primers <- filter(final.primers, orientation == "FORWARD") %>% 
  select(c(locus, SEQUENCE)) %>% arrange(locus) %>%
  mutate(loc.name = paste0(locus, ".R1"))
rev.primers <- filter(final.primers, orientation == "REVERSE") %>% 
  select(c(locus, SEQUENCE)) %>% arrange(locus)%>%
  mutate(loc.name = paste0(locus, ".R2"))
for (i in 1:nrow(fwd.primers)){
  fname <- paste0("data-raw/fasta/", description, ".", fwd.primers$loc.name[i], ".fasta")
  cat(paste0("> ", fwd.primers$loc.name[i], "\n", fwd.primers$SEQUENCE[i], "\n"), file = fname)
  fname <- paste0("data-raw/fasta/", description, ".", rev.primers$loc.name[i], ".fasta")
  cat(paste0("> ", rev.primers$loc.name[i], "\n", rev.primers$SEQUENCE[i], "\n"), file = fname)
}

# write amplicons to fasta
targets <- filter(targets, locus %in% final.primers$locus) %>% 
  left_join(filter(final.primers, orientation == "FORWARD")) %>%
  select(c(locus, fragment, start, PRODUCT_SIZE))
targets$start <- as.numeric(targets$start)
targets$PRODUCT_SIZE <- as.numeric(targets$PRODUCT_SIZE)
amplicons <- 
  mutate(targets, amp = substr(x = fragment, start = (start+1), stop = (start+PRODUCT_SIZE))) %>%
  select(c(locus, PRODUCT_SIZE, amp)) %>%
  mutate(amp.length = nchar(amp))
amplicon.fasta.file <- paste0("data-raw/fasta/", description, ".amplicons.fasta")
if(file.exists(amplicon.fasta.file)){
  print("fasta file already exists; delete it or change the name")
} else{
  for (i in 1:nrow(amplicons)){
    cat(paste0(">", amplicons$locus[i], "\n", amplicons$amp[i], "\n"), 
        file = amplicon.fasta.file, append = TRUE)
  }
}

# check to see if amplicons have extra Ns. If so, that means the part of the
# reference genome they mapped to had Ns, so the mapping might be unreliable
Ns.per.amp <- do.call(rbind, lapply(1:nrow(amplicons), function(i){
  Ns <- length(strsplit(amplicons$amp[i], split = "N")[[1]])-1
  return(data.frame(amplicons$locus[i], Ns))
})) %>% data.frame()
names(Ns.per.amp) <- c("locus", "Ns")
Ns.per.amp <- left_join(filtered.windows, Ns.per.amp)
Ns.per.amp$extra.Ns <- Ns.per.amp$Ns - Ns.per.amp$n.snps
length(which(Ns.per.amp$extra.Ns > 0))
# 141 amplicons have extra Ns; probably want to ditch those