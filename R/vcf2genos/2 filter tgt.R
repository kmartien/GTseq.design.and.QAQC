library(tidyverse)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(swfscMisc)
source("R/0 Load vcf2genos functions.R")
#source("R/filter.on.AR.and.minDP.r")
#source("R/Freebayes.miscalls.r")
#source("R/multiplot_pdf.r")

description = "GTseq.val.all.new"
date = format(Sys.time(), "%Y%b%d")

# Note that this assumes that the date on the saved file you're loading is today's
# date. If it was created previously, replace the file name below with the
# correct one.
load(paste0("data/", description, "_", date, "_vcf2tgt.rda"))

# filter out genotypes with "in between" allelic ratios (can't be called as het 
# or homo) and those that don't meet minimum depth (minDP) requirements specified
# in loc_ann. Also adds AR and minDP values as new columns to tgt. The column 
# 'locus' in loc_ann should be a factor
tgt.filtered <- filter.tgt(tgt, loc.ann)

########## CHECK FOR MISCALLS ##########
# Check for any freebayes miscalls (1-minAR) to 1 ratio heterozygotes or (1-maxAR) to maxAR
# ratio homozygotes. Use this to help correct leftover miscalls. Remember this list is not 
# comprehensive, you still need to review the plots for miscalls.

FB.miscalls <- freebayes.miscalls(tgt.filtered, out.file = paste("results/", description,"_", date, "_freebayes_miscalls.csv", sep = ""))

tgt$miscalled <- FALSE

for (i in 1:nrow(tgt)) {
  print(i)
  if (!is.na(tgt$GT[i])) {
    if (tgt$haplo.x[i] != tgt$haplo.y[i]) { #called as heterozygote
      if (tgt$ratio.x[i] >= (1-tgt$minAR[i])) { #should be REF allele homo
        tgt$miscalled[i] <- TRUE
        tgt$haplo.y[i] <- tgt$haplo.x[i]
      }
      if (tgt$ratio.y[i] >= (1-tgt$minAR[i])) { #should be ALT allele homo
        tgt$miscalled[i] <- TRUE
        tgt$haplo.x[i] <- tgt$haplo.y[i]
      }
    }
    if (tgt$haplo.x[i] == tgt$haplo.y[i] && !is.na(tgt$GT[i])) { #called as homozygote
      if ((1-tgt$maxAR[i]) >= tgt$ratio.x[i] & tgt$maxAR[i] <= tgt$ratio.x[i]) {
        tgt$miscalled[i] <- TRUE
        tgt$haplo.x[i] <- tgt$haplo.1[i]
        tgt$haplo.y[i] <- tgt$haplo.2[i]
      }
    }
  }
}

#next three lines for testing purposes
plot.limits <- data.frame(max = 500, mid = 100, min = 50)
plot.limits <- data.frame(max = sample(c(1000,500,750), length(unique(tgt$locus)), replace = TRUE), mid = 100, min = 50)
rownames(plot.limits) <- unique(tgt$locus)

plot.allele.ct.by.locus(tgt, locus_AR_minDP, fn.prefix = paste0("results/", description, "_", date, "_recalled"))

no.missing.genos.tgt <- tgt[-which(is.na(tgt$haplo.x)),]
missing.data.ind <- data.frame(table(no.missing.genos.tgt$Indiv)) %>%
  mutate(missing = (243-Freq)/243)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$missing > .1))

missing.data.ind <- separate_wider_delim(missing.data.ind, cols = "labID",
                                         delim = "_M", names = c("LabID","trash"))
missing.data.ind <- select(missing.data.ind, -trash)
data("gtseq.prod.summary")
gtseq.prod.summary$LabID <- paste0("z0", zero.pad(gtseq.prod.summary$LABID))

missing.data.ind <- left_join(missing.data.ind, gtseq.prod.summary)
pdf(file = "results/missing.data.pdf")
g <- ggplot(missing.data.ind, aes(x = On.Target.Reads, y = missing)) + geom_point()
g
dev.off()

missing.data.loc <- data.frame(table(no.missing.genos.tgt$locus)) %>%
  mutate(missing = 178-Freq)
length(which(missing.data.loc$missing > 17))

save(tgt, FB.miscalls, locus_AR_minDP, file = paste0("data/", description, "_", date, "_filtered-tgt.rda"))
