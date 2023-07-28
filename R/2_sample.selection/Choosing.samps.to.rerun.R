library(dplyr)
library(Mnov.GTseq.data)
library(swfscMisc)

data("gtseq.prod.summary")
data("gtseq.val.genos")
data("GTseq.samps.final")
data("id.key")
id.key$LABID <- paste0("z0", zero.pad(id.key$LABID))
gtseq.val.genos$LABID <- paste0("z0", zero.pad(gtseq.val.genos$LABID))

gtseq.smry.all <- rbind(prod.sum, gtseq.val.genos[,1:7])
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

conc <- read.csv("data-raw/sample.concentration.csv")
conc$LABID <- paste0("z0", zero.pad(conc$LABID))
conc.animals <- filter(id.key, LABID %in% conc$LABID ) %>%
  filter(id.type == "ANIMAL.id") %>% select(c(LABID, ANIMALID = alt.id)) %>%
  left_join(conc)
conc[-which(conc$LABID %in% conc.animals$LABID),]
names(conc.animals)[2] <- "AnimalID"

GTseq.samps.final <- left_join(GTseq.samps.final, conc.animals, by = "AnimalID")
