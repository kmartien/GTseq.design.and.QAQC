CCES.dat <- filter(MS54, USE.x == "CCES")
x.labs <- c("2018-07-04", "2018-07-17", "2018-08-02", "2018-08-10", "2018-09-08", 
            "2018-09-21","2018-10-22")

g <- ggplot(data = CCES.dat, aes(x = julian.date, y = gt.20.aligned)) + 
  geom_boxplot() + scale_x_discrete(labels = x.labs, breaks = format(as.Date(x.labs), "%j")) +
  labs(x = "Collection date", y = "Num. loci genotyped") +
  theme(text = element_text(size = 20))  

tiff("results-raw/CCES.success.vs.date.tiff", height = 900, width = 1800)
g
dev.off()
