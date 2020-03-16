setwd('D:Google Drive/Study/1. Semester an Uni-Tuebingen/Microarray Bioinformatik/HA6')
SingleExpression <- read.table(file = 'SingleExpressionCELfiles4A6.tsv', header = T, sep = "\t", nrows = 44041, 
                               row.names = 1) 

#scatter.smooth(LRMA1$A, LRMA1$M, col='red',degree = 1)
#scatter.smooth(LRMA1$A, LRMA1$M, col='red',degree = 2)

#LRMAnormalized <- LinearRegression(LRMA1$A, LRMA2$M)

par(mfcol = c(2,3))
lo1 = loess(LRMA1$M~LRMA1$A, span = 1/16, degree = 1)
lo2 = loess(LRMA2$M~LRMA2$A, span = 1/16, degree = 1)
lo3 = loess(LRMA3$M~LRMA3$A, span = 1/16, degree = 1)
lo4 = loess(LRMA4$M~LRMA4$A, span = 1/16, degree = 1)
lo5 = loess(LRMA5$M~LRMA5$A, span = 1/16, degree = 1)
lo6 = loess(LRMA6$M~LRMA6$A, span = 1/16, degree = 1)

fA1 <- predict(lo1, LRMA1$A)
fA2 <- predict(lo2, LRMA2$A)
fA3 <- predict(lo3, LRMA3$A)
fA4 <- predict(lo4, LRMA4$A)
fA5 <- predict(lo5, LRMA5$A)
fA6 <- predict(lo6, LRMA6$A)

M_normalized1 = LRMA1$M - fA1
M_normalized2 = LRMA2$M - fA2
M_normalized3 = LRMA3$M - fA3
M_normalized4 = LRMA4$M - fA4
M_normalized5 = LRMA5$M - fA5
M_normalized6 = LRMA6$M - fA6

plot(LRMA1$A, M_normalized1, pch='.')
plot(LRMA2$A, M_normalized2, pch='.')
plot(LRMA3$A, M_normalized3, pch='.')
plot(LRMA4$A, M_normalized4, pch='.')
plot(LRMA5$A, M_normalized5, pch='.')
plot(LRMA6$A, M_normalized6, pch='.')

# Judging from the produced plots through lowess and linear regression, it is clear that the plots produced through lowess
# show more equilibrium in the polarity of the group of points. Therefore, I claim that lowess is in this case better.

baselinerank <- matrix(rank(SingleExpression$Baseline.Array))
secondlinerank <- matrix(rank(SingleExpression$Second.Array))
ranklist <- cbind(baselinerank, secondlinerank)
SE <- cbind(SingleExpression, ranklist)
colnames(SE)[3:4] <- c("baserank", "secondrank")

k = nrow(SE)
tl = 0.003
th = 0.008

invagenebase = c()
invagenesecond = c()
size = rbind(NULL)
for (i in c(1:k)) {
  di = abs(SE[i, 3] - SE[i, 4]) / k
  rm = (SE[i, 3] + SE[i, 4]) / 2
  if(((di <= th) & (rm >= k/2)) | ((di <= tl) & (rm <= k/2))){
    invagenebase = rbind(invagenebase, SE[i, 1])
    invagenesecond = rbind(invagenesecond, SE[i, 2])
  }
  size = rbind(size, nrow(invagenebase))
}

invagene = cbind(invagenebase, invagenesecond)
colnames(invagene)[1:2] <- c("base", "second")

plot(c(1: nrow(size)), size, xlab = "Times of Iteration", ylab = "Size of the Set", type = "l", col = "red", 
     main = "Times of Iteration X Size of the Set")

invalr <- LinearRegression(invagene[,1], invagene[,2])
complr <- LinearRegression(SE$Baseline.Array, SE$Second.Array)

plot(SE$Baseline.Array, complr$modifiedchannel, pch=5, main = "Complete Normalized Data")
abline(a=complr$a, b=complr$b, col="blue", lty=3, lwd=3)

normalizedsecond <- as.matrix((SE$Second.Array - invalr$a) / invalr$b)

plot(SE$Baseline.Array, normalizedsecond[,1], pch=5, 
     main = "Normalized Data through linear regression of invariant genes")
abline(a=invalr$a, b=invalr$b, col="blue", lty=3, lwd=3)

# It is clear that the linear regression based on invariant genes would give us a better visualized outcome
# because of the equilibrium of expression values with lower fluctuation.