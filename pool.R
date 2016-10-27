# Combine allele frequencies from two samples (can be obtained with cleanCounts.R script)
# and calculate relative frequencies along chromosomes

library(zoo)
library(ggplot2)
library(grid)
library(scales)

##########################################
# import allele counts
##########################################

ctrlPool <- read.csv("control_counts.csv")
mutPool <- read.csv("mutant_counts.csv")


##########################################
# functions
##########################################

# remove positions with 0 counts for specific genotype in 2 datasets
zero <- function(x){
  x <- x[x[,3] + x[,5] != 0,]
  x <- x[x[,4] + x[,6] != 0,]
}

# pool neighboring snps; window size is set in number of snps
roll.chr <- function(x,window){
  for(i in 3:6){
    x <- cbind(x,rollapply(x[,i],
                           width=window,
                           FUN=sum,
                           by=window,
                           fill=NA))
  }
  names(x)[7:10] <- c(paste("sum.", names(x)[3:6], sep=""))
  x <- x[complete.cases(x), c(1, 2, 7:10)]
  return(x)
}

# calculate genotype frequencies
stats <- function(x){
  percCvi_B <- x[,6] / (x[,5] + x[,6])
  percCvi_A <- x[,4] / (x[,3] + x[,4])
  enrich <- (percCvi_B - percCvi_A) / percCvi_A
  cbind(x,
        percCvi_A = percCvi_A,
        percCvi_B = percCvi_B,
        enrich = enrich)
}

# smooth neighboring frequencies using the median
medianing <- function(x, window){
  mean.enrich <- rollapply(x$enrich, width = window, FUN = median, fill = NA)
  mean.ctrlPool <- rollapply(x$percCvi_A, width = window, FUN = median, fill = NA)
  mean.mutPool <- rollapply(x$percCvi_B, width = window, FUN = median, fill = NA)
  cbind(x,
        mean.enrich,
        mean.ctrlPool,
        mean.mutPool)
}

# arrange plots
vplayout <- function (x, y) {
  viewport(layout.pos.row = x, layout.pos.col = y)
}


##########################################
# merge datasets
##########################################

ds1 <- merge(ctrlPool[,c(1, 2, 7, 8)], mutPool[,c(1, 2, 7, 8)],
             by = c("CHR", "POS"), sort = FALSE,
             suffixes = c("ctrlPool", "mutPool"))

names(ds1)[1:2] <- c("Chr", "pos")

##########################################
# quality filter
##########################################

#remove NA's
ds1.cc <- ds1[complete.cases(ds1[,3:6]),]

#remove positions with no Ler or Cvi reads (probably incorrect genotype calls)
# (not advisable if the sample has a low coverage)
ds1.z <- zero(ds1.cc)

##########################################
# calculate frequencies along chromosomes 
##########################################

# split in each of the 5 A. thaliana chromosomes
data <- list()
for(i in 1:5){
  x <- ds1.z[substr(ds1.z$Chr,4,4) == i,]
  x <- x[order(x$pos),]
  data[[i]] <- x
}

# pool neighbouring snps using rolling window,
# calculate frequencies and average along neighbour region frequencies

for(i in 1:5){
  data[[i]] <- roll.chr(data[[i]], 50)
  data[[i]] <- stats(data[[i]])
  data[[i]] <- medianing(data[[i]], 100)
}

##########################################
# plot
##########################################

data.m <- rbind(data[[1]],
                data[[2]],
                data[[3]],
                data[[4]],
                data[[5]])

# use only subset of data to reduce image size
set.seed(42)
subdata <- data.m[sort(sample(1:dim(data.m)[1], size = 2000, replace = FALSE)),]
subdata$Chr <- factor(subdata$Chr,
                      labels=c("chr 1", "chr 2", "chr 3", "chr 4", "chr 5"))

# use Mbp for positions along chromosomes
subdata$pos <- subdata$pos / 1e6

# plot frequencies of Cvi (B reads) in each of the sets
a <- ggplot(data = subdata)+
  facet_grid(. ~ Chr, scales = "free_x", space = "free_x")+
  geom_ribbon(aes(x = pos, ymin = mean.ctrlPool, ymax = mean.mutPool), fill = "grey92")+
  geom_line(aes(pos, mean.ctrlPool, color = 'ctrlPool'))+
  geom_line(aes(pos, mean.mutPool, color = 'mutPool'))+
  scale_y_continuous(name = 'B allele frequency', labels = percent) +
  scale_x_continuous(labels = NULL,name = NULL) +
  scale_colour_manual(name = 'sets', values = c(ctrlPool = "blue", mutPool = "red"))+
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = 'top',
        strip.text.x = element_text(size = 12))

# plot relative enrichment of Cvi in mutant pool over the control pool
b <- ggplot(data = subdata, aes(pos, mean.enrich))+
  facet_grid(. ~ Chr, scales = "free_x", space = "free_x")+
  geom_hline(yintercept = 0, color = "grey") +
  geom_area(fill = 'gray98') +
  geom_line(color = 'grey30') +
  scale_y_continuous(name = 'B enrichment in mutant pool', labels = percent) +
  scale_x_continuous(name = 'Mbp') +
  theme_minimal()+
  theme(strip.text.x = element_blank())

# create grid to place plots
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1,)))

# print plots
print(a, vp = vplayout(1,1))
print(b, vp = vplayout(2,1))
