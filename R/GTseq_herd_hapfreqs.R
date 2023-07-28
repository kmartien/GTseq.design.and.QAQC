#I'VE UPDATED SOME STRATA IN THIS R INSTANCE, BUT NEED TO WRITE TO THE RDA AND
#UPDATE THE REPO. AS SOON AS GITHUB STARTS WORKING...

library(dplyr)
library(Mnov.GTseq.data)

data("AS.gtseq")
data("Mnov.strata")
data("ind.info")
data("GTseq.samps.final")

gtseq.ind.haps <- left_join(AS.gtseq, select(ind.info, c(LABID, SEX, HAP))) %>%
  left_join(select(Mnov.strata, c(LABID, wint.area, feed.area))) %>% 
  left_join(select(GTseq.samps.final,c(LABID, USE))) 

known.herd <- filter(gtseq.ind.haps, wint.area != "") %>% 
  filter(feed.area != "")%>%
  mutate(stratum = paste(wint.area,feed.area,sep="-"))

other.samps <- gtseq.ind.haps[-which(gtseq.ind.haps$LABID %in% known.herd$LABID),] %>%
  mutate(stratum = paste0(wint.area,feed.area))

all.df <- rbind(known.herd,other.samps) %>% select(c(LABID,stratum,HAP))

all.gtype <- df2gtypes(all.df, ploidy = 1, sequences = seqs)

all.successes.df <- left_join(all.df, read.sum) %>%
  filter(On.Target.Reads >= 20000)

