# Extracts target fragments from the chosen loci. Each fragment is 259bp long
# and is centered on the center of a SNP window. Generates an input file for
# Primer3 and calls Primer3 to generate primers. Results are formatted and 
# returned for further processing.

get.primer3.primers <- function(description, tidy.vcf, final.windows, 
                                chosen.loc.seqs, config.folder, fname){
  if(file.exists(fname)) {
    print(paste0("File ", fname, " already exists. Delete it, or choose a different file name."))
    return(NULL)
  }
  target.fragments <- do.call(rbind,lapply(names(chosen.loc.seqs), function(n){
    #get positions of all SNPs, including those I'm not targeting, and change them to n's
    snps <- filter(tidy.vcf, CHROM == n) %>% select(c(CHROM, POS, REF, ALT))
    snp.pos <- unique(snps$POS)
    n.4.snps <- as.character(chosen.loc.seqs[[which(names(chosen.loc.seqs)==n)]])
    for (i in 1:length(snp.pos)){
      ref.nuc <- tolower(snps[which(snps$POS==snp.pos[i])[1],3])
      if(n.4.snps[snp.pos[i]]==ref.nuc) n.4.snps[snp.pos[i]] <- 'n' else stop(paste("The nucleotide in chosen.loc.seqs does not match what's expected from the vcf.", n, i,sep=" "))
    }
    # calculate start, stop, and center positions of 259 bp fragments centered 
    # on my chosen windows
    target.windows <- filter(final.windows, CHROM == n) %>%
      select(start, stop) %>% mutate(length = stop-start+1) %>%
      mutate(center = (start+round((stop-start)/2))) %>%
      mutate(start.pos.in.frag = (start - center + 130)) %>%
      mutate(stop.pos.in.frag = (stop - center + 130))
    
    # create the input file for Primer3
    frags <- do.call(rbind, lapply(1:nrow(target.windows), function(w){
      s <- n.4.snps[(target.windows$center[w]-129):(target.windows$center[w]+129)]
      template <- paste(toupper(s), collapse = "")
      target.name <- paste0(n, "_", target.windows$start[w],"_", 
                            target.windows$stop[w])
      write(paste0("SEQUENCE_ID=", n, "_", target.windows$start[w],"_", 
                  target.windows$stop[w]),file=fname,append = TRUE)
      write(paste0("SEQUENCE_TEMPLATE=", template, collapse = ""), 
            file = fname, append = TRUE)
      write(paste0("SEQUENCE_TARGET=", target.windows$start.pos.in.frag[w],
                   ",", target.windows$length[w]), file = fname, append = TRUE)
      write("PRIMER_PRODUCT_SIZE_RANGE=90-143", file = fname, append = TRUE)
      write("PRIMER_MIN_GC=30", file = fname, append = TRUE)
      write("PRIMER_MAX_GC=70", file = fname, append = TRUE)
      write(paste0("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=", config.folder), file = fname, append = TRUE)
      write("=", file = fname, append = TRUE)
      return(c(target.name, template))
    }))
    return(frags)
  }))
  names(target.fragments) <- c("window", "target")
  primer3.res <- system2(command = "primer3_core",
                         args = c(fname),
                         stdout = TRUE, stderr = TRUE)
  primer.deets <- data.frame(do.call(rbind, strsplit(primer3.res, split = "=")))
  seq.id.lines <- c(which(primer.deets[,1] == "SEQUENCE_ID"), (nrow(primer.deets)+1))
  primer.deets$locus <- do.call(c, lapply(1:(length(seq.id.lines)-1), function(i){
    rep(primer.deets$X2[seq.id.lines[i]], times = seq.id.lines[i+1] - seq.id.lines[i])
  }))
  rows2keep <- c(grep("PRIMER_LEFT", primer.deets$X1), 
                 grep("PRIMER_RIGHT", primer.deets$X1), 
                 grep("PRIMER_PAIR", primer.deets$X1))
  primer.deets <- primer.deets[rows2keep,]
  primer.deets <- do.call(rbind, lapply(strsplit(primer.deets$X1, split = "_"), function(i){
    if(length(i) > 4) i[4] <- paste(i[4:length(i)], collapse = "_")
    if (length(i) < 4) i <- c(i, NA)
    return(i[1:4])
  })) %>% cbind(primer.deets)
  names(primer.deets) <- c("primer", "orientation", "index", "name","fullname","value", "locus")
  primer.deets$orientation[which(primer.deets$orientation %in% c("PAIR", "LEFT"))] <- "FORWARD"
  primer.deets$orientation[which(primer.deets$orientation == "RIGHT")] <- "REVERSE"
  primer.deets <- select(primer.deets, c(locus, orientation, index, name, value)) %>%
    filter(index != "NUM")
  start.pos.and.length <- filter(primer.deets, is.na(name))
  start.pos.and.length <- bind_cols(start.pos.and.length, do.call('rbind',strsplit(start.pos.and.length$value,",",fixed=TRUE))) %>%
    select(-c(name, value))
  names(start.pos.and.length)[4:5] <- c("start", "length")
  start.pos.and.length.long <- pivot_longer(start.pos.and.length, cols = c(start, length))
  primer.deets <- bind_rows(primer.deets, start.pos.and.length.long) %>%
    filter(!is.na(name)) %>% filter(!name == "PENALTY")
  primer.deets <- pivot_wider(primer.deets)
  return(list(targets=target.fragments, primer.deets=primer.deets))
}
