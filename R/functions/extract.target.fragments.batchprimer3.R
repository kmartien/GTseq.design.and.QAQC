# Extracts target fragments from the chosen loci. Each fragment is 259bp long
# and is centered on the center of a SNP window. Square brackets are added
# before the first SNP in the window and after the last one to identify for 
# Primer3 the portion of the sequence that should be amplified by the primers
# it designs.

extract.target.fragments <- function(tidy.vcf, final.windows, chosen.loc.seqs){
  target.frags <- do.call('c',lapply(names(chosen.loc.seqs), function(n){
    #get positions of all SNPs, including those I'm not targeting, and change them to n's
    snps <- filter(tidy.vcf, CHROM == n) %>% select(c(CHROM, POS, REF, ALT))
    snp.pos <- unique(snps$POS)
    n.4.snps <- as.character(chosen.loc.seqs[[which(names(chosen.loc.seqs)==n)]])
    for (i in 1:length(snp.pos)){
      ref.nuc <- tolower(snps[which(snps$POS==snp.pos[i])[1],3])
      if(n.4.snps[snp.pos[i]]==ref.nuc) n.4.snps[snp.pos[i]] <- 'n' else stop(paste("The nucleotide in chosen.loc.seqs does not match what's expected from the vcf.", n, i,sep=" "))
    }
    #extract 250 bp fragments centered on my chosen windows
    target.windows <- filter(final.windows, CHROM == n) %>%
      select(start, stop) %>% mutate(length = stop-start) %>%
      mutate(center = (start+round((stop-start)/2))) %>%
      mutate(start.pos.in.frag = (start - center + 130)) %>%
      mutate(stop.pos.in.frag = (stop - center + 130))
    # add brackets before the first SNP in the window and after the last one
    frags <- lapply(1:nrow(target.windows), function(w){
      s <- n.4.snps[(target.windows$center[w]-129):(target.windows$center[w]+129)]
      s2 <- paste(c(toupper(s[1:(target.windows$start.pos.in.frag[w] - 1)]), '[', 
                    toupper(s[target.windows$start.pos.in.frag[w]:target.windows$stop.pos.in.frag[w]]),
                    ']', toupper(s[(target.windows$stop.pos.in.frag[w]+1):259])))
    })
    names(frags) <- paste(n,target.windows$start, target.windows$stop, sep="_")
    return(frags)
  }))
  return(target.frags)
}