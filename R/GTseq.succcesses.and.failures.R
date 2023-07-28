library(dplyr)
library(Mnov.GTseq.data)

data("gtseq.prod.summary")
data("gtseq.val.genos")
data("GTseq.samps.final")
data("id.key")

gtseq.smry.all <- rbind(gtseq.prod.summary, gtseq.val.genos[,1:7])
gtseq.animals <- filter(id.key, LABID %in% gtseq.smry.all$LABID ) %>%
   filter(id.type == "ANIMAL.id") %>% select(c(LABID, ANIMALID = alt.id))
  
# Confirm missing AnimalIDs are all non-humpback
missing.samples <- GTseq.samps.final[-which(GTseq.samps.final$AnimalID %in% gtseq.animals$ANIMALID),]
names(missing.samples)[c(1,2)] <- c("ANIMALID","LABID")
if (length(which(missing.samples$GTseq.stratum %in% 
                 c("Blue", "brydes", "Fin", "Gray", "minke"))) < nrow(missing.samples)) 
{print("STOP!")} else {
    gtseq.animals <- rbind(gtseq.animals, missing.samples[,c(2,1)])
  }

read.sum <- do.call(rbind,lapply(gtseq.animals$LABID, function(i){
  inds <- filter(gtseq.smry.all, LABID == i)
  if (nrow(inds) == 1) {cbind(inds[1,c(1,3,4)], 1)} else {
    c(i, colSums(inds[,3:4]), nrow(inds))
  }
}))
names(read.sum)[4] <- "dupe"
read.sum$success <- read.sum$On.Target.Reads >= 20000

###   read.sum$On.Target.Reads SHOWS TOTAL READS PER SAMPLE, WITH READS FROM DUPLICATE SAMPLES COMBINED
save(read.sum, file = "results/GTseq.sample.success.summary.rda")

successes <- filter(read.sum, success == TRUE)
inds <- left_join(successes, GTseq.samps.final)
table(inds$GTseq.stratum)
