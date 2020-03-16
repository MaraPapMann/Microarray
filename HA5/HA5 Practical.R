# Reading Files
setwd("D:/Google Drive/Study/1. Semester an Uni-Tuebingen/Microarray Bioinformatik/HA5")
affy <- (read.table('affy_data.tsv', header = TRUE, row.names = 1))[1:54851,]

# Programming functions
LinearRegression = function(x, y){
  xm <- mean(x)
  ym <- mean(y)
  b <- sum((x-xm) * (y-ym)) / sum((x-xm)**2)
  a <- ym - b*xm
  return(list(a=a, b=b, gerade=function(z) a+b*z,
              modifiedchannel=(y-a)/b))
}

# Calculating linear regression of each pair probe
LR1 <- LinearRegression(affy$E13R021_P001.CEL, affy$E13R021_P002.CEL)
LR2 <- LinearRegression(affy$E13R021_P001.CEL, affy$E13R021_P005.CEL)
LR3 <- LinearRegression(affy$E13R021_P001.CEL, affy$E13R021_P011.CEL)
LR4 <- LinearRegression(affy$E13R021_P002.CEL, affy$E13R021_P005.CEL)
LR5 <- LinearRegression(affy$E13R021_P002.CEL, affy$E13R021_P011.CEL)
LR6 <- LinearRegression(affy$E13R021_P005.CEL, affy$E13R021_P011.CEL)

# Plotting original datas
par(mfcol=c(2,3))

plot(affy$E13R021_P001.CEL, affy$E13R021_P002.CEL, pch=20)
abline(a=LR1$a, b=LR1$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P001.CEL, affy$E13R021_P005.CEL, pch=20)
abline(a=LR2$a, b=LR2$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P001.CEL, affy$E13R021_P011.CEL, pch=20)
abline(a=LR3$a, b=LR3$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P002.CEL, affy$E13R021_P005.CEL, pch=20)
abline(a=LR4$a, b=LR4$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P002.CEL, affy$E13R021_P011.CEL, pch=20)
abline(a=LR5$a, b=LR5$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P005.CEL, affy$E13R021_P011.CEL, pch=20)
abline(a=LR6$a, b=LR6$b, col='blue', lwd=3, lty=1)

# Plotting normalized data

plot(affy$E13R021_P001.CEL, LR1$modifiedchannel, pch=20)
abline(a=LR1$a, b=LR1$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P001.CEL, LR2$modifiedchannel, pch=20)
abline(a=LR2$a, b=LR2$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P001.CEL, LR3$modifiedchannel, pch=20)
abline(a=LR3$a, b=LR3$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P002.CEL, LR4$modifiedchannel, pch=20)
abline(a=LR4$a, b=LR4$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P002.CEL, LR5$modifiedchannel, pch=20)
abline(a=LR5$a, b=LR5$b, col='blue', lwd=3, lty=1)

plot(affy$E13R021_P005.CEL, LR6$modifiedchannel, pch=20)
abline(a=LR6$a, b=LR6$b, col='blue', lwd=3, lty=1)

# Processing data with log2 function
affylog2 <- log2(affy)

# Loading library affy
library(affy)

# Programming MA calculating function
MA = function(x,y){
  M <- x-y
  A <- (x+y)/2
  M_mean <- mean(M)
  A_mean <- mean(A)
  b <- sum((M-M_mean) * (A-A_mean)) / sum((A-A_mean)**2)
  a <- M_mean - b*A_mean
  return(list(M=M, A=A, a=a, b=b, regression=function(z) a+b*z, modifiedchannel=(M-a)/b ))
}

# Calculating
LRMA1 <- MA(affylog2$E13R021_P001.CEL, affylog2$E13R021_P002.CEL)
LRMA2 <- MA(affylog2$E13R021_P001.CEL, affylog2$E13R021_P005.CEL)
LRMA3 <- MA(affylog2$E13R021_P001.CEL, affylog2$E13R021_P011.CEL)
LRMA4 <- MA(affylog2$E13R021_P002.CEL, affylog2$E13R021_P005.CEL)
LRMA5 <- MA(affylog2$E13R021_P002.CEL, affylog2$E13R021_P011.CEL)
LRMA6 <- MA(affylog2$E13R021_P005.CEL, affylog2$E13R021_P011.CEL)



# Plotting MA Plot with Linear Regression line
plot(LRMA1$A, LRMA1$M, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA1$a, b=LRMA1$b, col='red', lwd=1)

plot(LRMA2$A, LRMA2$M, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA2$a, b=LRMA2$b, col='red', lwd=1)

plot(LRMA3$A, LRMA3$M, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA3$a, b=LRMA3$b, col='red', lwd=1)

plot(LRMA4$A, LRMA4$M, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA4$a, b=LRMA4$b, col='red', lwd=1)

plot(LRMA5$A, LRMA5$M, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA5$a, b=LRMA5$b, col='red', lwd=1)

plot(LRMA6$A, LRMA6$M, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA6$a, b=LRMA6$b, col='red', lwd=1)



# Plotting Normalized MA Plot with Linear Regression line
plot(LRMA1$A, LRMA1$modifiedchannel, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA1$a, b=LRMA1$b, col='red', lwd=1)

plot(LRMA2$A, LRMA2$modifiedchannel, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA2$a, b=LRMA2$b, col='red', lwd=1)

plot(LRMA3$A, LRMA3$modifiedchannel, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA3$a, b=LRMA3$b, col='red', lwd=1)

plot(LRMA4$A, LRMA4$modifiedchannel, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA4$a, b=LRMA4$b, col='red', lwd=1)

plot(LRMA5$A, LRMA5$modifiedchannel, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA5$a, b=LRMA5$b, col='red', lwd=1)

plot(LRMA6$A, LRMA6$modifiedchannel, xlab = 'A', ylab = 'M', pch=20)
abline(a=LRMA6$a, b=LRMA6$b, col='red', lwd=1)

# After linear regression normalization, the coordinate was reset according to the results of linear regression.
# The vertical values were adjusted to a normalized level.
# Linear regression can be directly applied to scatterplot, but not to MA-plot.
# M'= M-f(A)

