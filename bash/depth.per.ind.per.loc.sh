for FILE in *.bam; do samtools depth -b /Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.data/data-raw/Mnov.targetSNP.bed $FILE > ${FILE%.bam}.coverage; done