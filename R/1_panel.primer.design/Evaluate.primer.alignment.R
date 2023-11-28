library(dplyr)
library(Mnov.GTseq.data)
source("./R/functions/extract.number.of.hits.R")

#data("Mnov.GTSEEK.panel")

sam.folder.path <- "data-raw/sam.files/Paired.primers.mapped.to.genome"

primer.matches <- extract.number.of.hits(sam.folder.path)

primer.matches$in.panel <- primer.matches$Locus %in% Mnov.GTSEEK.panel$short

GTSEEK.primer.test <- read.csv("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/SNPs/GTSEEK_validation_run/HB-Whale_PrimerTest.csv")
primer.matches <- left_join(primer.matches, GTSEEK.primer.test)

save(primer.matches, file = "data/primer.num.matches.to.genome.rda")
write.csv(primer.matches, file = "results-raw/primer.num.matches.to.genome.csv")


#evaluate amplicon mapping to reference
cov.folder.path <-"./JA_DATA/amplicon_coverage/"
locnames<-modbed$ampname
samplenames<-read.table("JA_DATA/amplicon_coverage/rad1sampnames.txt")
samplenames<-samplenames$V1

ampcoverage<-data.frame(matrix(NA,nrow=417,ncol=1))
names(ampcoverage)<-"locus"
ampcoverage$locus<-locnames

for (name in samplenames) {  
  infile <- paste0(cov.folder.path, name, ".cov")
  cov <-read.table(infile)
  if (nrow(cov) == 0) { 
  print(infile, "is empty, skipping to the next file.")
    next
  }
  locus.cov<-cov[,c(1,6)]
  names(locus.cov)<-c("locus",name)
  ampcoverage<- left_join(ampcoverage,locus.cov)
                       
                       
}
  
ampcoverage<- ampcoverage %>% mutate(total = rowMeans(pick(where(is.numeric))))

#now for depth:
ampdepth<-data.frame(matrix(NA,nrow=417,ncol=1))
names(ampdepth)<-"locus"
ampdepth$locus<-locnames

for (name in samplenames) {  
  infile <- paste0(cov.folder.path, name, ".cov")
  cov <-read.table(infile)
  if (nrow(cov) == 0) { 
    print(infile, "is empty, skipping to the next file.")
    next
  }
  locus.depth<-cov[,c(1,7)]
  names(locus.depth)<-c("locus",name)
  ampdepth<- left_join(ampdepth,locus.depth)
  
  
}

ampdepth<- ampdepth %>% mutate(total = rowMeans(pick(where(is.numeric))))
min(ampdepth$total) #3.79
max(ampdepth$total) #12266.44

#try removing samples where more than half of the loci have zero coverage
samplesums<-colSums(ampdepth[,-1])
min(samplesums)
mean(samplesums)
quantile(samplesums)
fail<-which((samplesums<100))
fail<-names(fail)
ampdepthnofail<-ampdepth[, -which(names(ampdepth) %in% fail)]
ampdepthnofail<- ampdepthnofail%>% mutate(total = rowMeans(pick(where(is.numeric))))

#this didn't really change calculations

multimappers<-read.table("JA_DATA/loci_to_remove.txt")
multimappers<-multimappers$V1
ampdepth2<-ampdepthnofail[-which(ampdepthnofail$locus %in% multimappers),]
write.csv(ampdepth2,"amplicondepth_allsamps.csv")

hist(ampdepth2$total[ampdepth2$total<100])
summary(ampdepth2$total)
sum(ampdepth2$total > 500) #22 loci
sum(ampdepth2$total > 1.5*210) #28 loci, mean is 210
#if I drop those 28 loci, we are down to 373 loci
sum(ampdepth2$total > 75) #46 loci, we would drop down to 355 loci
sum(ampdepth2$total > 2*15) #84 loci, 2*median we would drop down to 317 loci

#removing loci with greater than 1.5*mean (28 loci)
ampdepth3<-ampdepth2[which(ampdepth2$total <= 210*1.5),] #372 loci

#make bed file

 finalloci<-as.data.frame(map_chr(str_split(ampdepth3$locus,":"), 1))
finalloci$start <- (map_chr(str_split(map_chr(str_split(ampdepth3$locus,":"), 2),"-"),1))
   finalloci$stop<-(map_chr(str_split(map_chr(str_split(ampdepth3$locus,":"), 2),"-"),2))
   write.table(finalloci, "filteredloci11.26.23.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
   
multimapprimers<-read.table("./JA_DATA/multihitblast.txt")
multimappers<-multimapprimers %>% separate(V2, into=c("locus","read"), sep = "\\.R")
#look for duplicate loci indiating both primer pairs sare multimapping
multimappers2<-multimappers %>% group_by(locus) %>% filter(n()>1) #32 F/R pairs
multimappers2$locus <-gsub("\\.",":", multimappers2$locus) #match with previous dataframes for rbind
multimappers2<-multimappers2 %>% left_join(x= multimappers2, y= locs, by = join_by(locus == ampname))

multimapprimersamps<-multimappers2 %>% filter(!locus %in% multimappers)
finalloci11.27.23<-locs %>% filter(!ampname %in% multimappers) %>% filter(!ampname %in% multimappers2$locus)
finalloci11.27.23 %>% summarise(sum(no.snps))
sum(finalloci11.27.23$no.snps) #1067 snps