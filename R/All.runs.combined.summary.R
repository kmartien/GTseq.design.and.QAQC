library(dplyr)
library(swfscMisc)
library(Mnov.GTseq.data)

data("GTseq.samps.final")
data("id.key")
id.key$LABID <- paste0("z0", zero.pad(id.key$LABID))

gtseq.smry.all <- read.csv("results-raw/samps.in.each.run.csv")
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

gtseq.smry.all <- left_join(gtseq.smry.all, gtseq.animals)
names(gtseq.smry.all)[11] <- "AnimalID"
gtseq.smry.all <- left_join(GTseq.samps.final, gtseq.smry.all, by = "AnimalID")

min.genos <- 276

successes <- filter(gtseq.smry.all, genos.10readsMin >= min.genos)
table(successes$GTseq.stratum)
table(gtseq.smry.all$GTseq.stratum)
to.abandon <- filter(gtseq.smry.all, genos.10readsMin <50)
table(to.abandon$GTseq.stratum)
