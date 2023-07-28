library(dplyr)
library(Mnov.GTseq.data)
library(ggplot2)

data("gtseq.prod.summary")
data("gtseq.val.genos")
data("GTseq.samps.final")
data("id.key")

gtseq.smry.all <- rbind(gtseq.prod.summary, gtseq.val.genos[,1:7])
load("data/plate.info.rds")

gtseq.smry <- left_join(gtseq.smry.all, plate.info) %>% filter(!is.na(LABID))
gtseq.smry$plate[which(is.na(gtseq.smry$plate))] <- "p"

plate.sum <- data.frame(do.call(rbind, lapply(unique(gtseq.smry$plate), function(p){
  samps <- filter(gtseq.smry, plate == p)
  colMeans(samps[,3:7])
})))
plate.sum <- data.frame(cbind("plate" = unique(gtseq.smry$plate), plate.sum))
plate.sum <- mutate(plate.sum, Off.Target.Reads = Raw.Reads - On.Target.Reads)

to.plot <- pivot_longer(plate.sum, cols = c(On.Target.Reads, Off.Target.Reads))
g <- ggplot(to.plot, aes(fill = name, y = value, x = plate)) + 
  geom_bar(position = "stack", stat = "identity")
pdf(file = "results/reads.per.plate.pdf")
g
dev.off()
