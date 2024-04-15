rm(list=ls())
library(seqinr)
library(vcfR)
library(tidyr)
library(dplyr)
library(tidyverse)
library(here)

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


#Try reading in the 472 500 bp regions I've selected from dataset as fasta
chosen.loc.seqs<-read.fasta(file="./JA_DATA/primerregions2.fa", set.attributes = TRUE)
#chosen.loc.seqs<-do.call('rbind',strsplit(names(chosen.loc.seqs),":",fixed=TRUE, cols(c("chromo","length")))[,1])
#chosen.loc.seqs<- chosen.loc.seqs %>% separate_wider_delim(names(chosen.loc.seqs, delim = ":", names= c("chromo","length")))

#vcf2 <- read.vcfR("JA_DATA/phased_2019HIRAD.gt.vcf", convertNA = T)
#vcf2<- v[which(tidy_gt$chromo==n),c(1,2,3,4)]
#try with mafs file which only has snps, vcf has nonvariant sites also
vcf<- read.table("JA_DATA/2019RAD_snpsonly.mafs")
colnames(vcf)<-vcf[1,]
vcf<-vcf[-1,]
tidy_gt<-vcf

#tidy_gt2 <- vcfR2tidy(vcf2, 
                     single_frame = TRUE, 
                     info_fields = c("DP"), #"AS", "AD", "DP", "GQ", "AC", "AN", "PRO", "PAO", "AB", "DPRA", "ODDS", "MQM", "MQMR"
                     format_fields = c("GT", "GL", "AD", "RO", "QR", "AO", "QA")) #"GQ", "AC", "DP", "MIN_DP"=NA
#tidy_gt2$dat$CHROM <- do.call('rbind',strsplit(tidy_gt2$dat$CHROM,".",fixed=TRUE))[,1]

seq.frags <- do.call('c',lapply(names(chosen.loc.seqs), function(n){
  #get positions of all SNPs, including those I'm not targeting, and change them to n's
  tempn<-as.data.frame(chosen.loc.seqs[[n]])
  startpos<-as.numeric(map_chr(str_split(map_chr(str_split(names(chosen.loc.seqs[n]),":"), 2),"-"),1)) #the name of the n at the time to only have the positions, and also add the chromo to it)
  tempn$position<-seq(startpos,(startpos+500))
  chromn<-map_chr(str_split(names(chosen.loc.seqs[n]),":"), 1)
  snps <- tidy_gt %>% filter(tidy_gt$chromo==chromn, tidy_gt$position >= startpos,tidy_gt$position <= (startpos+500))
  snp.pos <- unique(snps$position)
  for (i in 1:length(snp.pos)){
    tempn$x[tempn$position==snp.pos[i]] <- '[N]' 
  }
  
  frags<-  as.data.frame(toupper(tempn$x))
  names(frags)<- n
  return(frags)
  
  }))


#Write target amplicons to files.  BatchPrimer3 only accepts 500 amplicons at a time, so
#the amplicons are written to multiple files in groups of 500

lapply(1:length(seq.frags), function(i){
  fname = paste("target_amplicons_",(floor(i/500)+1),".fasta",sep="")
  write(paste(">",names(seq.frags)[i],sep=" "),file=fname,append = TRUE)
  write(paste(seq.frags[[i]][1:501],sep="",collapse=""),file=fname,append=TRUE)
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

#if only reading in one file, changed to csv inexcel first
primer.deets <- read.csv("JA_DATA/target_amplicons_1.4.csv")
#batchprimer3 only generated one primer set per locus so I can't rang by importance or heterozygosity?

#calculate how many snps inside the insert region for each locus:

#first add in columns for the chromo start position for the locus of interest, and the true start and end positions of the amplicon
primer.deets$locuschrom<-map_chr(str_split(primer.deets$Seq.ID,":"), 1)
primer.deets$locusstart<-as.numeric(map_chr(str_split(map_chr(str_split(primer.deets$Seq.ID,":"), 2),"-"),1))
primer.deets$locusend <-as.numeric(map_chr(str_split(map_chr(str_split(primer.deets$Seq.ID,":"), 2),"-"),2))
primer.deets$ampstart <-primer.deets$locusstart + primer.deets$Start -1
primer.deets$ampend <- primer.deets$ampstart + primer.deets$Prod.Size -1 

#get number of ns (snps) in each amplicon
Fonly<-  primer.deets %>% filter(Orientation == "FORWARD")
amplicon.snps<- as.data.frame(matrix(NA,(nrow(Fonly)),3))
colnames(amplicon.snps)<- c("locus","no.snps","rankorder")
#tidy_gt$position <-as.numeric(tidy_gt$position)

 for (i in 1:nrow(Fonly)) {
  
  amplicon.snps[i,1]<- Fonly$Seq.ID[i]
  amplicon.snps[i,2]<-tidy_gt %>% filter(tidy_gt$chromo == Fonly$locuschrom[i],tidy_gt$position >= Fonly$ampstart[i], tidy_gt$position <= Fonly$ampend[i]) %>% count()
  amplicon.snps[i,3]<-Fonly$Count[i]
}
  
#append snp count to primer.deets
primer.deets <- left_join(x=primer.deets, y=amplicon.snps, join_by("Seq.ID" == "locus" ,"Count" == "rankorder"), keep = FALSE)

primer.deets %>% filter(Orientation == "FORWARD") %>% filter(no.snps >0) %>% count()
primer.deets %>% filter(Orientation == "FORWARD") %>% filter(no.snps >0) %>% count(no.snps)
primer.deets %>% filter(Orientation == "FORWARD") %>% filter(no.snps >0) %>% summarise(sum(no.snps))

#211 loci with at least 1 snp, total of 342 SNPs, not amazing but better than what I've been getting? No snps in primer regions, all between 90-143bp product
#no.snps   n
#1       1 128
#2       2  55
#3       3  20
#4       4   5
#5       8   3


#ran batchprimer 3 again with the same data and conditions and it generated a different set (target_amplicons_1.2.csv)
#279 loci with at least 1 snp, total of 561 SNPs
#no.snps   n
#1       1 117
#2       2  85
#3       3  55
#4       4  13
#5       5   3
#6       6   3
#7       7   1
#8       8   1
#9       9   1



#running again with sam data BUT chaning product length to be between 170-220 (target_amplicons_1.3.csv)
#359 loci with at least 1 snp, total of 877 snps
#  no.snps   n
#1       1 102
#2       2  93
#3       3 104
#4       4  40
#5       5  11
#6       6   5
#7       7   1
#8       8   2
#9       9   1

#Other values to consider/
#the value for "any" is the calvulated score of the tendenc of a primer to bind to itself, scores ANY binding occurring within the entire primer sequence. the max is 8
#there is another score for tendeny of left primer to bind to righ tprimer with a max value of 8 also. Similar for the 3' binding (which would result in primer dimers and the default val is 3)
#See no reason to be concerned with the values I have currently for their tendency to form primer dimers, so moving on

#try setting primer 3 to output 5 primer sets per locus instead of 1 to sort that way. Find that even with multiple primer sets per locus the number of snps isn't that different
#find loci where none of the primer sets result in snps in insert:

primer.deets %>% group_by(Seq.ID) %>% summarise(max_snps = max(no.snps)) %>% filter(max_snps == 0) %>% count()
#52 loci have no primer sets that result in snps , leaving 420 loci

no_snp_loci<-primer.deets %>% group_by(Seq.ID) %>% summarise(max_snps = max(no.snps)) %>% filter(max_snps == 0) 

onlysnps_primerdeets<-primer.deets %>% filter(!Seq.ID %in% no_snp_loci$Seq.ID)

checkloci<-onlysnps_primerdeets %>% group_by(Seq.ID) %>% summarise(has_zero = any(no.snps == 0)) %>% filter(has_zero == TRUE)
checkprimers<-onlysnps_primerdeets %>% filter(Seq.ID %in% checkloci$Seq.ID) %>% filter(no.snps != 0)


#want to choose the primer pair for each locus that results in the highest number of snps, and if there are multiple with the same snp count choose the lowest count value ( highest rank)

unique_seq_ids <- unique(onlysnps_primerdeets$Seq.ID)
selected_rows <-data.frame()

for (i in 1:length(unique_seq_ids)) {
  print(paste0("working with seq_id", unique_seq_ids[i])) 
  subset_data <- onlysnps_primerdeets[onlysnps_primerdeets$Seq.ID == unique_seq_ids[i], ]  # Subset for each Seq.ID
  max_snps <- max(subset_data$no.snps)  # Find maximum no.snps within subset
  filtered_subset <- subset_data[subset_data$no.snps == max_snps, ]  # Subset where no.snps is highest

  if (nrow(filtered_subset) > 2) {
    min_count <- min(filtered_subset$Count)  # Find minimum count within subset
    final_subset <- filtered_subset[filtered_subset$Count == min_count, ]  # Subset where count is lowest
    selected_rows <- rbind(selected_rows, final_subset)  # Append to selected_rows
  } else {
    selected_rows <- rbind(selected_rows, filtered_subset)  # Append to selected_rows
  }
  
}

selected_rows %>% filter(Orientation == "FORWARD") %>% filter(no.snps >0) %>% count(no.snps)
#420 loci, snp distrbution: 1186 total snps
#no.snps   n
#1        1  70
#2        2 108
#3        3 144
#4        4  60
#5        5  21
#6        6   8
#7        7   3
#8        8   4
#9        9   1
#10      13   1

#note that SUPER_14:4263211-4263380, SUPER_17:21988648-21988821, SUPER_17:23218405-23218574, SUPER_19:7365224-7365724 ended up in here twice, removed from fasta file downstream


##SNPs that don't have another target SNP within 100,000 bp are good.  For those that do 
#have another target SNP within that distance, I need to check to decide which one to keep

locs <- selected_rows  %>% filter(Orientation == "FORWARD")
locs <-locs %>% group_by(locuschrom) %>% mutate(ampdist = ampstart - lag(ampend))
locs <- locs %>% mutate(CHECK = ifelse(ampdist <=10000, "CHECK","PASS"))
locs %>% filter(CHECK == "CHECK") %>% count()
locs$ampname <- paste0(locs$locuschrom,":",locs$ampstart,"-",locs$ampend)
#filter out strange duplicates that got in there?
locs <- locs %>% filter(!duplicated(ampname))
selected_rows <- selected_rows %>% filter(Seq.ID %in% locs$Seq.ID)
# There are 79 amplicons that are within 100000 bases of each other, should filter out? That woul dleave 341, or if we do 10k 


# Make bed file of amplicon start and end for extracting fasta sequence and blasting amplicons to each other
write.table(locs[,c("locuschrom","ampstart","ampend","Seq.ID")],"distfiltloci.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table((selected_rows %>% filter(Orientation == "FORWARD"))[,c("locuschrom","ampstart","ampend","Seq.ID")],"nodistfiltloci.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

modbed<-selected_rows %>% filter(Orientation == "FORWARD") %>% select("locuschrom","ampstart","ampend","Prod.Size") %>% mutate(ampname = paste0(locuschrom,":",ampstart,"-",ampend)) %>% mutate(amp0 = 0) %>% mutate(ampstoppos = 0 + Prod.Size)

write.table((selected_rows %>% filter(Orientation == "FORWARD"))[,c("locuschrom","ampstart","ampend","Prod.Size")],"nodistfiltloci.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(modbed[,5:7],"bedforgatk.bed",sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)



######################## Karen's filters ###############

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
