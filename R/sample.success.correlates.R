library(Mnov.GTseq.data)
library(strataG)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(psych)
library(viridis)
library(rfPermute)
library(randomForest)

data("AS.gtseq")
data("ind.info")
data("Mnov.strata")
data("gtseq.prod.summary")
data("gtseq.val.genos")
data("tissue.archive.data")
data("GTseq.samps.final")
data("id.key")

gtseq.smry.all <- rbind(gtseq.prod.summary, gtseq.val.genos[,1:7])
gtseq.animals <- filter(id.key, LABID %in% gtseq.smry.all$LABID ) %>%
  filter(id.type == "ANIMAL.id") %>% select(c(LABID, ANIMALID = alt.id)) %>%
  left_join(select(ind.info, c("LABID","HAP","SEX")), by = "LABID") %>%
  left_join(select(tissue.archive.data, c("LABID","YR","MO","DA","State","Country","Latitude","Longitude")), by = "LABID") %>%
  left_join(gtseq.smry.all, by = "LABID") %>%
  left_join(read.csv("data-raw/sample.concentration.csv")) %>%
  mutate(successful = ifelse(On.Target.Reads >= 20000,1,0)) %>%
  mutate(Dried.down = ifelse(Dried.down == "X", 1, 0)) %>%
  filter(Concentration <= 200)

# Plot all variables against each other and display correlation coefficients
discrete_palette <- viridis(n=4)[4:1]
title <- "Correlation plots"
pdf(file = "results/correlation.plots.Conc.lt.200.pdf")
pairs.panels(select(gtseq.animals, c(On.Target.Reads, Raw.Reads, Pct.On.Target, successful, Concentration, YR, Dried.down)), ellipses = FALSE, stars = TRUE,
             gap = 0,
             bg = "gray",
#             bg = discrete_palette,
             main = title, pch=21)
dev.off()

gtseq.animals$Raw.Reads <- as.numeric(gtseq.animals$Raw.Reads)
gtseq.animals$On.Target.Reads <- as.numeric(gtseq.animals$On.Target.Reads)
dups <- filter(gtseq.animals, LABID %in% names(which(table(gtseq.animals$LABID) > 1)))
ordered.dups <- dups[with(dups, order(LABID, Raw.Reads)),]
dups <- left_join(select(dups[seq(from= 1, to = 166, by = 2),], LABID, Raw.Reads, On.Target.Reads, Pct.On.Target),select(dups[seq(from = 2, to = 166, by = 2),], LABID, Raw.Reads, On.Target.Reads, Pct.On.Target), by = "LABID")
title <- "Duplicates correlation plots"
pdf(file = "results/dups.correlation.plots.Conc.lt.200.pdf")
pairs.panels(dups[,-1], ellipses = FALSE, stars = TRUE,
             gap = 0,
             bg = "gray",
             #             bg = discrete_palette,
             main = title, pch=21)
dev.off()

# Random forest analysis to predict 'success' based on other variables

success.df <- select(gtseq.animals, c("Raw.Reads", "Pct.On.Target", 
  "successful", "Concentration", "YR", "Dried.down")) %>%
  filter(!is.na(Concentration)) %>% filter(!is.na(YR))
success.df$successful <- factor(c("Failure", "Success")[success.df$successful + 1])
samp.size <- balancedSampsize(success.df$successful)
rf <- rfPermute(successful ~ ., success.df, sampsize = samp.size, replace = FALSE, ntree = 3000, nrep = 1000, importance = TRUE)
pdf("results/Variable importance.Conc.lt.200.pdf")
plotImportance(rf)
dev.off()

bonus.samps <- read.csv("data-raw/bonus.samples.csv") %>% left_join(gtseq.animals) %>%
  filter(!is.na(Raw.Reads))

# Still need to correct the collection info for the 4 samples where I had the wrong labID
other.sp.smry <- gtseq.smry.all[-which(gtseq.smry.all$LABID %in% gtseq.animals$LABID),] %>%
  left_join(GTseq.samps.final, by = "LABID") %>% 
  select(-c(AnimalID,Duplicate1, Duplicate2, Duplicate3, Herd, HAP, mito.haps, USE, Optimization.sample, Main.run)) 
