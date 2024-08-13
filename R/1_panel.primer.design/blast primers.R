library(tidyverse)

# 1) Download and install BLAST+ executable from here: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

# 2) Make sure your amplicons are in a single FASTA file in the current working directory

amplicon.fasta.file <- 'Bphy_amplicons.fasta'


# 3) Create a local blast database from this FASTA file

if(!dir.exists('amplicons_db')) dir.create('amplicons_db')
system2(
  "makeblastdb",
  paste(
    paste0("-in ", amplicon.fasta.file), 
    "-parse_seqids",
    "-dbtype nucl",
    "-out amplicons_db/amplicons_db"
  ),
  wait = TRUE
)


# 4) BLAST the same FASTA file against the database

blast.cols <- c(
  "qseqid", "sacc",
  "evalue", "bitscore", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "slen"
)
amplicon.blast <- system2(
  command = "blastn",  
  args = c(
    "-db amplicons_db/amplicons_db",
    paste0("-query ", amplicon.fasta.file),
    paste0("-outfmt '6 ", paste(blast.cols, collapse = " "), "'"), 
    "-task 'megablast'"
  ),
  wait = TRUE,
  stdout = TRUE
) %>%
  as_tibble() %>% 
  separate(col = value, into = blast.cols, sep = "\t", convert = TRUE) %>% 
  mutate(
    pct.length = 100 * length / slen,
    pct.mismatch = 100 * mismatch / length
  ) %>% 
  arrange(desc(pident), evalue)
  

# 5) Remove self matches and hits that are at less than 60% of the sequence length
#   (other filtering criteria could be added here...)

to.delete <- amplicon.blast %>% 
  filter(qseqid != sacc & pct.length >= 60) %>% 
  select(qseqid, sacc) |> 
  unlist() |> 
  unique()