library(seqinr)
library(tidyr)
library(dplyr)

GTseq.loci <- read.csv("data-raw/rename_chrom.csv", header = FALSE)
names(GTseq.loci) <- c("locus", "new.name")
GTseq.loci <- separate_wider_delim(GTseq.loci, cols = locus, delim = "_", names = c("CHROM", "POS"), cols_remove = FALSE)
GTseq.loci$POS <- as.numeric(GTseq.loci$POS)

# First read in the contigs from the reference genome that contain target loci then skip to line 42 
load(file="results-R/GTseq.loc.seqs.full.CHROM.Rdata")

###############################################################################
# If the above file doesn't exist, read in the reference genome, which I split 
# into five files, and select the contigs that contain target loci.  That way
# I don't have to keep the whole huge reference in memory
ref.genome <- read.fasta(file="/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.1.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.1 <- ref.genome[which(names(ref.genome) %in% GTseq.loci$CHROM)]
rm(ref.genome)
ref.genome <- read.fasta(file="/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.2.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.2 <- ref.genome[which(names(ref.genome) %in% GTseq.loci$CHROM)]
rm(ref.genome)
ref.genome <- read.fasta(file="/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.3.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.3 <- ref.genome[which(names(ref.genome) %in% GTseq.loci$CHROM)]
rm(ref.genome)
ref.genome <- read.fasta(file="/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.4.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.4 <- ref.genome[which(names(ref.genome) %in% GTseq.loci$CHROM)]
rm(ref.genome)
ref.genome <- read.fasta(file="/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.5.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.5 <- ref.genome[which(names(ref.genome) %in% GTseq.loci$CHROM)]
rm(ref.genome)

GTseq.loc.seqs <- c(chosen.loc.seqs.1,chosen.loc.seqs.2,chosen.loc.seqs.3,chosen.loc.seqs.4,chosen.loc.seqs.5)
save(GTseq.loc.seqs,file="data/GTseq.loc.seqs.full.CHROM.Rdata")
################################################################################

seq.frags <- do.call('c',lapply(names(GTseq.loc.seqs), function(n){
  seq <- as.character(GTseq.loc.seqs[[which(names(GTseq.loc.seqs)==n)]])
  target.snps <- GTseq.loci$POS[which(GTseq.loci$CHROM == n)]
  frags <- lapply(target.snps, function(p){
    seq[(p-129):(p+129)]
  })
  names(frags) <- paste(n,target.snps,sep="_")
  return(frags)
}))

#Write target amplicons to files
fname = paste("data-raw/fasta/Mnov_target_amplicons-old_names-all.fasta",sep="")
lapply(1:length(seq.frags), function(i){
  write(paste(">",names(seq.frags)[i],sep=" "),file=fname,append = TRUE)
  write(do.call(paste0, lapply(1:259, function(l){toupper(seq.frags[[i]][l])})), file = fname, append = TRUE)
})

new.names <- data.frame(locus = names(seq.frags)) %>% left_join(GTseq.loci)
names(seq.frags) <- new.names$new.name

fname = paste("data-raw/fasta/Mnov.target.amplicons.all.fasta",sep="")
lapply(1:length(seq.frags), function(i){
  write(paste(">",names(seq.frags)[i],sep=" "),file=fname,append = TRUE)
  write(do.call(paste0, lapply(1:259, function(l){toupper(seq.frags[[i]][l])})), file = fname, append = TRUE)
})

