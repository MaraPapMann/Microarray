setwd('C:/Users/isoch/Google ‘∆∂À”≤≈Ã/Study/1. Semester an Uni-Tuebingen/Microarray Bioinformatik/HA4')
ad <- read.table(file=file.choose(), row.names = 1, header = T, sep = "\t", stringsAsFactors = T)
affydata <- ad[1:54851,]

# Load Libraries
library(ggplot2)
library(gplots)

# Combine plots
par(mfcol=c(2,3))

# Produce pairwise scatterplots of the primary expression values of each pair of arrays
plot(affydata$E13R021_P001.CEL, affydata$E13R021_P002.CEL, pch=20)
abline(a=2500, b=1, col='blue',lwd=5, lty=3)
abline(a=-2500, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydata$E13R021_P001.CEL, affydata$E13R021_P005.CEL, pch=20)
abline(a=2500, b=1, col='blue',lwd=5, lty=3)
abline(a=-2500, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydata$E13R021_P001.CEL, affydata$E13R021_P011.CEL, pch=20)
abline(a=2500, b=1, col='blue',lwd=5, lty=3)
abline(a=-2500, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydata$E13R021_P002.CEL, affydata$E13R021_P005.CEL, pch=20)
abline(a=2500, b=1, col='blue',lwd=5, lty=3)
abline(a=-2500, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydata$E13R021_P002.CEL, affydata$E13R021_P011.CEL, pch=20)
abline(a=2500, b=1, col='blue',lwd=5, lty=3)
abline(a=-2500, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydata$E13R021_P005.CEL, affydata$E13R021_P011.CEL, pch=20)
abline(a=2500, b=1, col='blue',lwd=5, lty=3)
abline(a=-2500, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

# Logtransform
affydatalog2 <- log2(affydata)


# Produce scatterplots of the log-transformed expression values of each pair of arrays.
plot(affydatalog2$E13R021_P001.CEL, affydatalog2$E13R021_P002.CEL, pch=20)
abline(a=1, b=1, col='blue',lwd=5, lty=3)
abline(a=-1, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydatalog2$E13R021_P001.CEL, affydatalog2$E13R021_P005.CEL, pch=20)
abline(a=1, b=1, col='blue',lwd=5, lty=3)
abline(a=-1, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydatalog2$E13R021_P001.CEL, affydatalog2$E13R021_P011.CEL, pch=20)
abline(a=1, b=1, col='blue',lwd=5, lty=3)
abline(a=-1, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydatalog2$E13R021_P002.CEL, affydatalog2$E13R021_P005.CEL, pch=20)
abline(a=1, b=1, col='blue',lwd=5, lty=3)
abline(a=-1, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydatalog2$E13R021_P002.CEL, affydatalog2$E13R021_P011.CEL, pch=20)
abline(a=1, b=1, col='blue',lwd=5, lty=3)
abline(a=-1, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

plot(affydatalog2$E13R021_P005.CEL, affydatalog2$E13R021_P011.CEL, pch=20)
abline(a=1, b=1, col='blue',lwd=5, lty=3)
abline(a=-1, b=1, col='blue',lwd=5, lty=3)
abline(a=0, b=1, col='red',lwd=3, lty=2)

# Produce MA Plots
library(affy)

ma.plot(rowMeans(affydatalog2[,c(1,2)]), affydatalog2[,1] - affydatalog2[,2], cex = 1)
ma.plot(rowMeans(affydatalog2[,c(1,3)]), affydatalog2[,1] - affydatalog2[,3], cex = 1)
ma.plot(rowMeans(affydatalog2[,c(1,4)]), affydatalog2[,1] - affydatalog2[,4], cex = 1)
ma.plot(rowMeans(affydatalog2[,c(2,3)]), affydatalog2[,2] - affydatalog2[,3], cex = 1)
ma.plot(rowMeans(affydatalog2[,c(2,4)]), affydatalog2[,2] - affydatalog2[,4], cex = 1)
ma.plot(rowMeans(affydatalog2[,c(3,4)]), affydatalog2[,3] - affydatalog2[,4], cex = 1)

# It is obvious that from a MA plot we can see the difference between the expression level of both probes more directly and
# clearly. From the Median and IQR of each MA plot we may tell the expression level generally in every single experiment.
# In conclusion, MA plot should be more appropriate to visualize the data. 

# Remove all the plots
dev.off()
