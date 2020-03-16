
#  a)

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("stringr", repos='http://cran.us.r-project.org')

source("http://bioconductor.org/biocLite.R")
biocLite("HELP")
library("HELP")

setwd('C:/Users/isoch/Google ‘∆∂À”≤≈Ã/Study/1. Semester an Uni-Tuebingen/Microarray Bioinformatik/HA2')
library(Biostrings)
library(stringr)

#  b)

#Data reading
o1 = readDNAStringSet("oligos2017.fasta")
o2 <- data.frame(o1)

#  c)

#Function programming
MT <- function(x){
  numa = sum(str_count(x, "A"))
  numt = sum(str_count(x, "T"))
  numc = sum(str_count(x, "C"))
  numg = sum(str_count(x, "G"))
  Tm = 2*(numa + numt) + 4*(numc + numg)
  return(Tm)
}

Tm <- apply(o2[,1], 1, MT)

print('The minimum of Tm is')
print(min(Tm))
print('The maximum of Tm is')
print(max(Tm))
print('The average of Tm is')
print(mean(Tm))
print('The median of Tm is')
print(median(Tm))

boxplot(Tm,col='green',
        pars=list(outcol='red'),
        main='Melting Temperature',
        xlab='Oligo',
        ylab='Temperature/Celcius degree')

#According to the boxplot, these oligos are apparently not a good choice for a microarray 
#because the melting temperature of the oligos are way too different, which will greatly influence(negatively) 
#the sensitivity of the microarray.

#  d)
# Random Oligo Generating

library(stringi)
bases <- c("A", "T", "C", "G")
randomOligos <- function(n){
  do.call(paste0, replicate(n, sample(bases, 1000, TRUE), FALSE))
}


AverageTmRandom <- function(x){
  data1 = data.frame(randomOligos(x))
  TmArray = apply(data1, 1, MT)
  AverageTm = mean(TmArray)
  return(AverageTm)
}

array1 = (1:300)
array2 = data.frame(array1)
AverageTm300 = apply(array2, 1, AverageTmRandom)

plot(AverageTm300, main='Average of Melting Temperature of Randomly Generated Oligos', xlab='Length of Oligo / Bp', 
     ylab='Average of Melting Temperature / Celsius Degree')

# According to this plot, I suggest that a length between 50 and 100 should provide a proper temperature for pizza baking...