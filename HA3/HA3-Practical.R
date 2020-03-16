# Read Data
setwd("C:/Users/isoch/Google ‘∆∂À”≤≈Ã/Study/1. Semester an Uni-Tuebingen/Microarray Bioinformatik/HA3")
Pe <- as.data.frame(read.table("probeExp.tsv", header = TRUE, row.names = NULL))

# Calculate Median
medianPe <- aggregate(Pe[,-1], list(Pe[,1]), median)

# Calculate Variance of each gene across all experiments
varPe <- apply(medianPe[,-1], 1, var)

# Making results into a matrix
Mvar <- matrix(varPe, nrow = 499, ncol = 1, byrow = TRUE)

VaR <- cbind(medianPe[,1], Mvar)

write.table(VaR, file='Variance.tsv', quote=FALSE, sep='\t', col.names = NA)
