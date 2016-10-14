# Retrieve publicly available snp data for the Cvi-0 and Ler-1 accessions of Arabidopsis thaliana,
# merge and outputs a reformated .csv file

################################
#import snp files
################################

Ler <- read.table('http://1001genomes.org/data/MPI/MPISchneeberger2011/releases/current//Ler-1/Marker/Ler-1.SNPs.TAIR9.txt')
# more info can be found at http://1001genomes.org/data/MPI/MPISchneeberger2011/releases/current//README

Cvi <- read.table('http://signal.salk.edu/atg1001/data/Salk/quality_variant_filtered_Cvi_0.txt')
# more info can be found at http://signal.salk.edu/atg1001/data/README.txt

Ler_pos <- paste("chr", Ler[,2], "_", Ler[,3], sep = "")
Cvi_pos <- paste(Cvi[,2], Cvi[,3], sep = "_")

################################
#create snp table
################################

snp <- rbind(cbind(pos = Ler_pos,
                   Ler[, c(1, 4, 5)]),
             cbind(pos = Cvi_pos,
                   Cvi[, c(1, 4, 5)]))
snp$pos <- as.character(snp$pos)

#remove common snps
snp.m <- snp[!snp$pos %in% snp$pos[duplicated(snp$pos)],]

################################
# export to file
################################

write.csv(snp.m, file = "snpm.csv")
