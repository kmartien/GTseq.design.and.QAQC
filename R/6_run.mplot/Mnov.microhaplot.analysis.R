library(vcfR)
library(microhaplot)

run.label <- "RunMS51"
vcf <- "final.recode"

# for your dataset: customize the following paths
 sam.path <- "data-raw/sam.files/RunMS51"
 label.path <- "data-raw/mplot_labels/RunMS51.label.txt"
 vcf.path <- "vcf/RunMS51.targetSNPs.recode.vcf"
# app.path <- "~/Shiny/microhaplot"

#sam.path <- paste0("data-raw/sam.files/", run.label)
#label.path <- file.path("data-raw/mplot_labels/", paste0(run.label, ".label.txt"))
#vcf.path <- paste0("vcf/", run.label, ".targetSNPs.recode.vcf")
out.path <- "results-R/microhaplot"
app.path <- "/Users/Shared/KKMDocuments/Documents/Github.Repos/Shiny/microhaplot"

vcf <- read.vcfR(vcf.path)
#locus.ignore.file <- read.csv(paste0("microhaplot/",run.label, ".locus_annotation.csv"))

# I've prepped the data, so can just jump straight to running the Shiny app
haplo.read.tbl <- prepHaplotFiles(run.label = paste0(run.label,".allsamps"),
                                  sam.path = sam.path,
                                  out.path = out.path,
                                  label.path = label.path,
                                  vcf.path = vcf.path,
                                  app.path = app.path,
                                  n.jobs = 2) # assume running on dual core

runShinyHaplot(app.path)
