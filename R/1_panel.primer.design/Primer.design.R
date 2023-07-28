rm(list=ls())
library(seqinr)
library(vcfR)
library(tidyr)
library(dplyr)
library(tidyverse)


setwd("~/Documents/Karen/Structure/Humpbacks/SNPs/Mnov_ddRAD_noFlash")
load("Selected.SNPs.top10pct.Rdata")

#First read in the reference genome, which I split into five files, and select the contigs that contain target loci.  That way
#I don't have to keep the whole huge reference in memory

###  If you've already done this part, you can skip over it by reading in the file saved on line ~40 (save(chosen.loc.seqs,file="chosen.loc.seqs.Rdata"))
ref.genome <- read.fasta(file="~/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.1.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.1 <- ref.genome[which(names(ref.genome) %in% loc.sum$locus)]
save(chosen.loc.seqs.1,file="chosen.loc.seqs.1.Rdata")
rm(ref.genome)
ref.genome <- read.fasta(file="~/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.2.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.2 <- ref.genome[which(names(ref.genome) %in% loc.sum$locus)]
save(chosen.loc.seqs.2,file="chosen.loc.seqs.2.Rdata")
rm(ref.genome)
ref.genome <- read.fasta(file="~/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.3.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.3 <- ref.genome[which(names(ref.genome) %in% loc.sum$locus)]
save(chosen.loc.seqs.3,file="chosen.loc.seqs.3.Rdata")
rm(ref.genome)
ref.genome <- read.fasta(file="~/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.4.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.4 <- ref.genome[which(names(ref.genome) %in% loc.sum$locus)]
save(chosen.loc.seqs.4,file="chosen.loc.seqs.4.Rdata")
rm(ref.genome)
ref.genome <- read.fasta(file="~/Documents/Karen/Structure/Humpbacks/data/Reference.genomes/GCA_004329385.1.5.fna")
names(ref.genome) <- do.call('rbind',strsplit(names(ref.genome),".",fixed=TRUE))[,1]
chosen.loc.seqs.5 <- ref.genome[which(names(ref.genome) %in% loc.sum$locus)]
save(chosen.loc.seqs.5,file="chosen.loc.seqs.5.Rdata")
rm(ref.genome)

chosen.loc.seqs <- c(chosen.loc.seqs.1,chosen.loc.seqs.2,chosen.loc.seqs.3,chosen.loc.seqs.4,chosen.loc.seqs.5)
save(chosen.loc.seqs,file="chosen.loc.seqs.Rdata")

#Read in the filtered VCF file that resulted from the SNP discovery phase.  Change all SNPs to Ns, select the sequence flanking 
#each target locus, and write it to a file formatted to be read by BatchPrimer3

vcf <- read.vcfR("vcf/Mnov_filtered_SNPs.vcf", convertNA = T)

tidy_gt <- vcfR2tidy(vcf, 
                     single_frame = TRUE, 
                     info_fields = c("DP"), #"AS", "AD", "DP", "GQ", "AC", "AN", "PRO", "PAO", "AB", "DPRA", "ODDS", "MQM", "MQMR"
                     format_fields = c("GT", "GL", "AD", "RO", "QR", "AO", "QA")) #"GQ", "AC", "DP", "MIN_DP"=NA
tidy_gt$dat$CHROM <- do.call('rbind',strsplit(tidy_gt$dat$CHROM,".",fixed=TRUE))[,1]

seq.frags <- do.call('c',lapply(names(chosen.loc.seqs), function(n){
  #get positions of all SNPs, including those I'm not targeting, and change them to n's
  snps <- tidy_gt$dat[which(tidy_gt$dat$CHROM==n),c(1,2,4,5)]
  snp.pos <- unique(snps$POS)
  n.4.snps <- as.character(chosen.loc.seqs[[which(names(chosen.loc.seqs)==n)]])
  for (i in 1:length(snp.pos)){
    ref.nuc <- tolower(snps[which(snps$POS==snp.pos[i])[1],3])
    if(n.4.snps[snp.pos[i]]==ref.nuc) n.4.snps[snp.pos[i]] <- 'n' else stop(paste("The nucleotide in chosen.loc.seqs does not match what's expected from the vcf.", n, i,sep=" "))
  }
  #extract 129 bps on either side of each of my chosen snps
  target.snps <- chosen.locs$pos[which(chosen.locs$locus == n)]
  frags <- lapply(target.snps, function(p){
    n.4.snps[(p-129):(p+129)]
  })
  names(frags) <- paste(n,target.snps,sep="_")
  return(frags)
}))

#Write target amplicons to files.  BatchPrimer3 only accepts 500 amplicons at a time, so
#the amplicons are written to multiple files in groups of 500
lapply(1:length(seq.frags), function(i){
  fname = paste("Mnov_target_amplicons_",(floor(i/500)+1),".fasta",sep="")
  write(paste(">",names(seq.frags)[i],sep=" "),file=fname,append = TRUE)
  write(paste(c(toupper(seq.frags[[i]][1:129]), '[N]', toupper(seq.frags[[i]][131:259])),sep="",collapse=""),file=fname,append=TRUE)
  write("\n",file=fname,append=TRUE)
  names(seq.frags)[i]
})

###########################################################################
#
#         Upload files one at a time into BatchPrimer3 
# (https://probes.pw.usda.gov/cgi-bin/batchprimer3/batchprimer3.cgi)
#       Run with defaults except Product length from 90 to 143, GC content 
#       between 25% and 75%.  Download full results as .zip file, extract
#       zip into project directory, and rename each results folder as
#       "Mnov_primer_results_#". Save each summary file (end in ..tmp.txt) 
#       as a .csv file
#
###########################################################################

#Read in the details of the primers designed by BatchPrimer3.  They're in multiple
#files due to 500 locus limit on BatchPrimer3
primer.deets <- do.call('rbind',lapply(1:12, function(i){
  fname <- list.files(path = paste("Mnov_primer_results_",i,sep=""),pattern="tmp.csv")
  tab <- read.csv(file=paste("Mnov_primer_results_",i,"/",fname,sep=""),as.is = TRUE)
}))
primer.deets[,1:2] <- do.call("rbind",strsplit(primer.deets$Seq.ID,split="_",fixed=TRUE))
names(primer.deets)[1:2] <- c("locus","pos")
primer.deets$pos <- as.numeric(primer.deets$pos)

save(primer.deets, file="primer.deets.Rdata")

#I'm hoping to end up with 500 loci, so start by selecting the top ranked 350 from each the 
#RF importance ranking and the heterozygosity ranking.  The join of those two sets is my working list
top.imp <- subset(imp, sum.rank >= imp$sum.rank[dim(imp)[1]-350], select=c("locus","sum.rank"))
top.het <- subset(loc.rank, rank.sum >= loc.rank$rank.sum[dim(loc.rank)[1]-350])
top.locs.string <- data.frame(unique(c(top.imp$locus,top.het$locus)))
names(top.locs.string) <- "locus"
top.locs.string <-  left_join(top.locs.string,top.imp)
top.locs.string <- left_join(top.locs.string, top.het)
nrow(top.locs.string)
top.locs <- data.frame(do.call("rbind",strsplit(top.locs.string$locus,split="_")),stringsAsFactors = FALSE)
top.locs <- top.locs[-2]
names(top.locs) <- c("locus","pos")
top.locs$pos <- as.numeric(top.locs$pos)
top.locs <- cbind(top.locs,top.locs.string[-1])

#count the number of primer pairs successfully designed for each locus.  Remove loci that don't have primers
top.locs$n.primer.pairs <- as.numeric(sapply(1:nrow(top.locs), function(i){
  primer.deets %>% filter(locus == top.locs$locus[i]) %>% filter(pos == top.locs$pos[i]) %>% dplyr::summarise(n=dplyr::n()/2)
}))
top.locs <- top.locs %>% filter(n.primer.pairs > 0)

top.locs <- top.locs %>% 
  group_by(locus) %>%
  arrange(pos)

loc.sum <- top.locs %>% 
  group_by(locus) %>% 
  dplyr::summarise(
    n = dplyr::n(),
    max.pos = max(pos),
    min.pos = min(pos)
  ) %>% 
  mutate(length = max.pos - min.pos)

#SNPs that don't have another target SNP within 100,000 bp are good.  For those that do 
#have another target SNP within that distance, I need to check to decide which one to keep
loc.lists <- lapply(unique(top.locs$locus), function(x){
  locs <- top.locs %>% filter(locus==x) 
  if (nrow(locs) == 1) return(list(locs.to.keep=locs,locs.to.check=NULL)) else {
    locs$min.dist <- sapply(1:nrow(locs), function(i) {
      min(abs((locs$pos-locs$pos[i])[-i]))
    })
    locs.to.keep <- locs[which(locs$min.dist >=100000),]
    locs.to.check <- locs[which(locs$min.dist < 100000),]
    return(list(locs.to.keep=locs.to.keep,locs.to.check=locs.to.check))
  }
})
names(loc.lists) <- unique(top.locs$locus)

l.2.keep <- do.call('rbind', lapply(loc.lists, function(l){l$locs.to.keep}))
l.2.check <- do.call('rbind', lapply(loc.lists, function(l){l$locs.to.check}))
write.csv(l.2.check,file="locs.2.check.csv")

####################################################################################
#
#  I went through locs.2.check.csv manually and decided which position(s) to keep for each
#  contig such that all chosen loci are at least 100,000 bp apart.  Where there were multiple
#  target positions within 100bp of each, I tried to select primers that captured all of the
#  target SNPs.  There were 9 cases where no primers captured them all and I tried to redesign
#  primers by bracketing a larger target (including all target SNPs).  Only 3 of those 
#  succeeded in producing primers.  In a few cases I noted that a primer other than 
#  the first one needed to be kept in order to capture all of the SNPs or because of
#  problems with the first primer pair
#
####################################################################################

# read back in checked loci with only those using still included in the .csv file
checked.locs <- read.csv(file="Checked.locs.csv",as.is = TRUE)
locs.to.keep <- rbind(l.2.keep,l.2.check[checked.locs[,1],])

redesigned.primer.deets <- read.csv(file="Mnov_redesigned_primers/Mnov_redesigned_primers.csv",as.is = TRUE)
redesigned.primer.deets[,1:2] <- do.call("rbind",strsplit(redesigned.primer.deets$Seq.ID,split="_",fixed=TRUE))
names(redesigned.primer.deets)[1:2] <- c("locus","pos")
redesigned.primer.deets$pos <- as.numeric(redesigned.primer.deets$pos)

#default assumption is that we'll keep the first primer pair, since BatchPrimer3 thinks that's the best
locs.to.keep$pair.to.keep <- 1

primer.to.change <- checked.locs %>% filter(!is.na(pair.to.keep))

#Change which primer pair to keep for those loci that I manually chose a different one
for (i in 1:nrow(primer.to.change)) {
  p <- which(locs.to.keep$locus==primer.to.change$locus[i] & locs.to.keep$pos==primer.to.change$pos[i])
  locs.to.keep$pair.to.keep[p] <- primer.to.change$pair.to.keep[i]
}

#Unless I specified a primer pair to keep, keep the lowest-numbered one (BatchPrimer3 ranks them by quality) that 
#produces a product less than 130 bp
save(locs.to.keep, file="locs.to.keep.csv")
primers.to.order <- lapply(1:nrow(locs.to.keep),function(i){
  if(locs.to.keep$pair.to.keep[i] == 0) {  # 0 indicates a redesigned primer, so get its details from redesigned.primer.deets
    primer.pair <- redesigned.primer.deets %>% filter(locus == locs.to.keep$locus[i]) %>% filter(pos == locs.to.keep$pos[i]) %>% filter(Count == 1)
  } else {
    primer.pair <- primer.deets %>% filter(locus == locs.to.keep$locus[i]) %>% filter(pos == locs.to.keep$pos[i]) %>% filter(Count == Count[min(which(Prod.Size<130))])
  }
  return(list(fwd=primer.pair[1,],rev=primer.pair[2,]))
})

primers.to.order <- list(
  fwd = do.call('rbind', lapply(primers.to.order, function(l){l$fwd})),
  rev = do.call('rbind', lapply(primers.to.order, function(l){l$rev})))
to.keep <- which(!is.na(primers.to.order$fwd$locus))
primers.to.order$fwd <- slice(primers.to.order$fwd,to.keep)
primers.to.order$rev <- slice(primers.to.order$rev,to.keep)

plot.primer.specs(primers.to.order)

#I only ended up with 479 loci, so I need to select more.  Get the top 450 ranks from RF importance
#and heterozygosity, and remove the ones that were in the first round
next.imp <- subset(imp, sum.rank >= imp$sum.rank[dim(imp)[1]-450], select=c("locus","sum.rank"))
next.het <- subset(loc.rank, rank.sum >= loc.rank$rank.sum[dim(loc.rank)[1]-450])
next.locs.string <- data.frame(unique(c(next.imp$locus,next.het$locus)),stringsAsFactors = FALSE)
names(next.locs.string) <- "locus"
#this next line removes the loci that were in the first round
next.locs.string <- slice(next.locs.string, c(1:nrow(next.locs.string))[-which(next.locs.string$locus %in% top.locs.string$locus)])
next.locs.string <-  left_join(next.locs.string,next.imp)
next.locs.string <- left_join(next.locs.string,next.het)
nrow(next.locs.string)
next.locs <- data.frame(do.call("rbind",strsplit(next.locs.string$locus,split="_")),stringsAsFactors = FALSE)
next.locs <- next.locs[-2]
names(next.locs) <- c("locus","pos")
next.locs$pos <- as.numeric(next.locs$pos)
next.locs <- cbind(next.locs,next.locs.string[-1])

next.locs$n.primer.pairs <- as.numeric(sapply(1:nrow(next.locs), function(i){
  primer.deets %>% filter(locus == next.locs$locus[i]) %>% filter(pos == next.locs$pos[i]) %>% dplyr::summarise(n=dplyr::n()/2)
}))
next.locs <- next.locs %>% filter(n.primer.pairs > 0)

next.locs$old.new <- "new"
top.locs$old.new <- "old"
all.locs <- bind_rows(top.locs,next.locs)

#need to consider both rounds of loci when making sure my final list doesn't include loci within
#100,000 bp of each other.  Loci are labeled as 'old' (from the first round of top 350) and 
#'new' (from the next 100 top ranks) so that I always keep the old locus if there is a comflict
next.loc.lists <- lapply(unique(all.locs$locus), function(x){
  locs <- all.locs %>% filter(locus==x) 
  locs$min.dist <- NA
  if (nrow(locs) == 1) return(list(locs.to.keep=locs,locs.to.check=NULL)) else {
    locs$min.dist <- sapply(1:nrow(locs), function(i) {
      min(abs((locs$pos-locs$pos[i])[-i]))
    })
    locs.to.keep <- locs[which(locs$min.dist >=100000),]
    locs.to.check <- locs[which(locs$min.dist < 100000),]
    return(list(locs.to.keep=locs.to.keep,locs.to.check=locs.to.check))
  }
})
names(next.loc.lists) <- unique(all.locs$locus)

next.l.2.keep <- do.call('rbind', lapply(next.loc.lists, function(l){l$locs.to.keep}))
next.l.2.keep <- filter(next.l.2.keep,old.new == "new")

next.l.2.check <- do.call('rbind', lapply(next.loc.lists, function(l){l$locs.to.check}))
write.csv(next.l.2.check,file="next.locs.2.check.csv")

# read back in checked loci with only those using still included in the .csv file
next.checked.locs <- read.csv(file="Next.checked.locs.csv",as.is = TRUE)
next.locs.to.keep <- rbind(next.l.2.keep,next.l.2.check[next.checked.locs[,1],])

next.primers.to.order <- lapply(1:nrow(next.locs.to.keep),function(i){
  primer.pair <- primer.deets %>% filter(locus == next.locs.to.keep$locus[i]) %>% filter(pos == next.locs.to.keep$pos[i]) %>% filter(Count == Count[min(which(Prod.Size<130))])
  return(list(fwd=primer.pair[1,],rev=primer.pair[2,]))
})

next.primers.to.order <- list(
  fwd = do.call('rbind', lapply(next.primers.to.order, function(l){l$fwd})),
  rev = do.call('rbind', lapply(next.primers.to.order, function(l){l$rev})))
to.keep <- which(!is.na(next.primers.to.order$fwd$locus))
next.primers.to.order$fwd <- slice(next.primers.to.order$fwd,to.keep)
next.primers.to.order$rev <- slice(next.primers.to.order$rev,to.keep)

plot.primer.specs(next.primers.to.order, description = "Second.round")
write.csv(primers.to.order$fwd, file= "fwd.primers.to.order.csv")
write.csv(primers.to.order$rev, file= "rev.primers.to.order.csv")
write.csv(next.primers.to.order$fwd, file= "next.fwd.primers.to.order.csv")
write.csv(next.primers.to.order$rev, file= "next.rev.primers.to.order.csv")

save(primers.to.order,next.primers.to.order,file="primers.to.order.Rdata")

primer.seqs <- rbind(cbind(paste(primers.to.order$fwd$locus,primers.to.order$fwd$pos,"F",sep="_"),primers.to.order$fwd$Seq),
                     cbind(paste(primers.to.order$rev$locus,primers.to.order$rev$pos,"R",sep="_"),primers.to.order$rev$Seq))
colnames(primer.seqs) <- c("primer.names","Seqs")

fname <- "primer.seqs.fasta"
lapply(1:nrow(primer.seqs), function(x){
  write(paste(">",primer.seqs[x,1]),file=fname,append = TRUE)
  write(primer.seqs[x,2], file=fname, append=TRUE)
})
