library(Mnov.GTseq.data)
library(swfscMisc)
library(dplyr)

data("AS.gtseq")
data("ind.info")
data("Mnov.strata")
data("gtseq.prod.summary")
data("gtseq.val.genos")
data("tissue.archive.data")
data("GTseq.samps.final")
data("id.key")
id.key$LABID <- paste0("z0", zero.pad(id.key$LABID))
gtseq.val.genos$LABID <- paste0("z0", zero.pad(gtseq.val.genos$LABID))
ind.info$LABID <- paste0("z0", zero.pad(ind.info$LABID))
tissue.archive.data$LABID <- paste0("z0", zero.pad(tissue.archive.data$LABID))

conc <- read.csv("data-raw/sample.concentration.csv")
conc$LABID <- paste0("z0", zero.pad(conc$LABID))

gtseq.smry.all <- rbind(prod.sum, gtseq.val.genos[,1:7])
gtseq.animals <- filter(id.key, LABID %in% gtseq.smry.all$LABID ) %>%
  filter(id.type == "ANIMAL.id") %>% select(c(LABID, ANIMALID = alt.id)) %>%
  left_join(select(ind.info, c("LABID","HAP","SEX")), by = "LABID") %>%
  left_join(select(tissue.archive.data, c("LABID","YR","MO","DA","State","Country","Latitude","Longitude")), by = "LABID") %>%
  left_join(gtseq.smry.all, by = "LABID") %>%
  left_join(conc) %>%
  mutate(successful = ifelse(On.Target.Reads >= 20000,1,0)) %>%
  mutate(Dried.down = ifelse(Dried.down == "X", 1, 0)) 

gtseq.animals$OptOrProd <- "Prod"
gtseq.animals$OptOrProd[which(gtseq.animals$Sample %in% gtseq.val.genos$Sample)] <- "Opt"

#test.samps <- read.csv("data-raw/Mnov GTseq test samples.csv")
write.csv(gtseq.animals, file = "results-raw/choosing.test.samples.csv")
