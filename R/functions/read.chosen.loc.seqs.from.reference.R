library(seqinr)
library(vcfR)
library(tidyr)
library(dplyr)
library(tidyverse)

read.chosen.loc.seqs.from.reference <- function(final.windows, ref.genome.location, chosen.locs.file = "data/chosen.locs.seqs.rda"){

  # I split my reference genome into five files since the whole thing is so large.
  # I read in each file, select the contigs that contain target loci, then remove 
  # that file and move on to the next.  That way I don't have to keep the whole 
  # huge reference in memory
  
  ref.genome <- read.fasta(file= paste0(ref.genome.location, "1.fna"))
  names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
  chosen.loc.seqs.1 <- ref.genome[which(names(ref.genome) %in% final.windows$CHROM)]
  save(chosen.loc.seqs.1,file="chosen.loc.seqs.1.Rdata")
  rm(ref.genome)
  ref.genome <- read.fasta(file= paste0(ref.genome.location, "2.fna"))
  names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
  chosen.loc.seqs.2 <- ref.genome[which(names(ref.genome) %in% final.windows$CHROM)]
  save(chosen.loc.seqs.2,file="chosen.loc.seqs.2.Rdata")
  rm(ref.genome)
  ref.genome <- read.fasta(file= paste0(ref.genome.location, "3.fna"))
  names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
  chosen.loc.seqs.3 <- ref.genome[which(names(ref.genome) %in% final.windows$CHROM)]
  save(chosen.loc.seqs.3,file="chosen.loc.seqs.3.Rdata")
  rm(ref.genome)
  ref.genome <- read.fasta(file= paste0(ref.genome.location, "4.fna"))
  names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
  chosen.loc.seqs.4 <- ref.genome[which(names(ref.genome) %in% final.windows$CHROM)]
  save(chosen.loc.seqs.4,file="chosen.loc.seqs.4.Rdata")
  rm(ref.genome)
  ref.genome <- read.fasta(file= paste0(ref.genome.location, "5.fna"))
  names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
  chosen.loc.seqs.5 <- ref.genome[which(names(ref.genome) %in% final.windows$CHROM)]
  save(chosen.loc.seqs.5,file="chosen.loc.seqs.5.Rdata")
  rm(ref.genome)
  
  chosen.loc.seqs <- c(chosen.loc.seqs.1,chosen.loc.seqs.2,chosen.loc.seqs.3,chosen.loc.seqs.4,chosen.loc.seqs.5)
  save(chosen.loc.seqs,file="data/chosen.loc.seqs.rda")
  return(chosen.loc.seqs)
}