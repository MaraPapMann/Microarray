
library(httr)
library(magrittr)
library(plotly, quietly = T)

# Read Data
setwd("C:/Users/isoch/Google ‘∆∂À”≤≈Ã/Study/1. Semester an Uni-Tuebingen/Microarray Bioinformatik/HA3")
AQ <- read.delim("ArrayQuality.tsv")

# RLE Calculate
splitted <- split(AQ, AQ$Identifier)

RLE <- sapply(names(splitted), function(n) {
	dt <- splitted[[n]][,-1]
	med <- apply(dt, 1, median)
	res <- dt - med
}, simplify = F)

IQR <- sapply(names(splitted), function(n) {
	dt <- RLE[[n]]
	quantiles <- apply(dt, 2, quantile)
	iqr <- quantiles[4,] - quantiles[2,]
}) %>% t

MED <- sapply(names(splitted), function(n) {
	dt <- RLE[[n]]
	med <- apply(dt, 2, median)
}) %>% t


# plot
MED = data.frame(MED)
MED.stack <- stack(MED)

pt <- ggplot(MED.stack) + 
	geom_boxplot(aes(x=ind, y=values)) +
	labs(x="Array", y="Median of RLE") + 
	theme_bw()
ggsave("ArrayQuality.png", pt)


# NUSE Calculate
NUSE <- sapply(names(splitted), function(n) {
  dt <- splitted[[n]][,-1]
  sd <- apply(dt, 2, sd)
  len <- apply(dt, 2, length)
  se <- matrix(sd / len, nrow = 1, ncol = 5)
  medse <- apply(se, 1, median)
  nuse <- se / medse
}, simplify = T)

rownames(NUSE) = c("A1", "A2", "A3", "A4", "A5")
# Sorting Matrix
t_NUSE <- t(NUSE)
NUSEframe <- data.frame(t_NUSE)
NUSE.stack <- stack(NUSEframe)

# Plot
pt2 <- ggplot(NUSE.stack) + 
  geom_boxplot(aes(x=ind, y=values)) +
  labs(x="Array", y="Median of NUSE") + 
  theme_bw()
ggsave("ArrayQuality2.png", pt2)
