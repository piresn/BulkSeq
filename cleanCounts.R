# Merge info from the file snpm.txt (known SNPs matrix) with the determined allele frequencies (from an out.vcf file)
# and outputs a csv file with allele counts

################################
#functions
################################

import <- function(x){
  imp <- read.table(x)[, c(1, 2, 4, 5, 6, 8)]
  names(imp) <- c("chrom", "coord", "ref.seq", "alt.seq", "qual", "info")
  imp$chrom <- tolower(imp$chrom)
  imp$info <- as.character(imp$info)
  
  # get DP4 info where counts for different alleles are stored
  # See http://samtools.sourceforge.net/mpileup.shtml
  a <- sapply(strsplit(as.character(imp$info), split = "DP4="), "[", 2)
  b <- sapply(strsplit(a, split=";"), "[", 1)
  c1 <- as.numeric(sapply(strsplit(b, split=","), "[", 1))
  c2 <- as.numeric(sapply(strsplit(b, split=","), "[", 2))
  c3 <- as.numeric(sapply(strsplit(b, split=","), "[", 3))
  c4 <- as.numeric(sapply(strsplit(b, split=","), "[", 4))
  
  #combine
  ref <- c1 + c2
  alt <- c3 + c4
  pos <- paste(imp[, 1], imp[, 2], sep="_")
  out <- cbind(pos, imp[, -6],ref, alt)
  out$pos <- as.character(out$pos)
  out$chrom <- factor(out$chrom)
  return(out)
}


snp.info <- function(x, snp.m = snp.m){
  origin <- snp.m$genot[match(x$pos, snp.m$pos, nomatch = NA)]
  ref.origin <- snp.m$ref[match(x$pos, snp.m$pos, nomatch = NA)]
  alt.origin <- snp.m$alt[match(x$pos, snp.m$pos, nomatch = NA)]
  bind <- cbind(x,
                origin = origin,
                ref.origin = ref.origin,
                alt.origin = alt.origin)
  out <- bind[complete.cases(bind),]
  return(out)
}


filter <- function(x, cut){
  cat('Initial number of SNPs:', dim(x)[1], '\n\n')
  
  #use only snps for which reference is ATGC
  filt1 <- x[x$ref.seq %in% c("A", "T", "G", "C"),]
  cat(dim(x)[1] - dim(filt1)[1], 'reference sequence values were not ATGC\n')
  
  #use only snps that match between sequenced and imported snps
  filt2 <- filt1[as.character(filt1$alt.seq) == as.character(filt1$alt.origin),]
  cat(dim(filt1)[1] - dim(filt2)[1], 'snps did not match between reference and vcf files\n')
  
  # remove high coverage based on expected values from poisson distribution
  cov <- filt2$ref + filt2$alt
  filt3 <- filt2[p.adjust(ppois(cov, lambda = median(cov), lower.tail = F),
                          method = "bonferroni") > cut,]
  temp <- dim(filt2)[1] - dim(filt3)[1]
  cat(temp, 'snps (', round(temp*100 / dim(filt2)[1], 2), '%) were above the expected',
      cut*100, '% percentile of a theoretical Poisson distribution\n\n')
  
  cat('Returned number of SNPs:', dim(filt3)[1])
  return(filt3)
}


format <- function(x){
  count <- with(x,
                data.frame(CHR = chrom,
                           POS = coord,
                           COV = ref + alt,
                           TAIR10_BASE = ref.seq,
                           LER_BASE = factor('N', levels = c('A', 'T', 'G', 'C')),
                           CVI_BASE = factor('N', levels = c('A', 'T', 'G', 'C'))
                ))
  count$LER_COUNT[x$origin == "Ler-1"] <- x$alt[x$origin == "Ler-1"]
  count$CVI_COUNT[x$origin == "Ler-1"] <- x$ref[x$origin == "Ler-1"]
  count$CVI_COUNT[x$origin == "Cvi_0"] <- x$alt[x$origin == "Cvi_0"]
  count$LER_COUNT[x$origin == "Cvi_0"] <- x$ref[x$origin == "Cvi_0"]
  count$LER_BASE[x$origin == "Ler-1"] <- x$alt.seq[x$origin == "Ler-1"]
  count$CVI_BASE[x$origin == "Cvi_0"] <- x$alt.seq[x$origin == "Cvi_0"]
  
  cat("Median coverage:", median(count$COV))
  cat("\n")
  cat("Median Cvi read frequency:", with(count, median(CVI_COUNT / COV)))
  cat("\n")
  cat("Median Ler read frequency:", with(count, median(LER_COUNT / COV)))
  cat("\n")
  cat("Total snps:", dim(count)[1])
  return(count)
}

################################
# run
################################

# import snps
snp.m <- read.table('snpm.txt', header = TRUE)

# import allele frequencies from vcf file
vcf <- import("out.vcf")

# combine vcf with known snp info
vcf.info <- snp.info(vcf, snp.m)

# filter out ambiguous and top-covered snps
# the cut parameter removes the reads which are above the top  0.1% percentil of a theoretical distribution
vcf.filter <- filter(vcf.info, cut = 0.001)

# reformat
counts <- format(vcf.filter)

################################
# export
################################

write.csv(counts, 'counts.csv', row.names = FALSE, quote = FALSE)
