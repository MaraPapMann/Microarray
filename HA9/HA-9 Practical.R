# Data Reading
library(gplots)
library(ComplexHeatmap)
library(cluster)

setwd("C:/Users/isoch/Google �ƶ�Ӳ��/Study/1. Semester an Uni-Tuebingen/Microarray Bioinformatik/HA8/material")

Spell_Man_Cell_Cycle_Fullset <- read.table("Spellman_cellcycle_alpha_completeexpressionmatrix.txt", header = T)

Name_list_Cell_Cycle_Gene <- read.table("Spellman_cellcycle_alpha_subsetcellcycle_geneIDs.txt")

S_cycle <- subset(Spell_Man_Cell_Cycle_Fullset, row.names(Spell_Man_Cell_Cycle_Fullset) %in% Name_list_Cell_Cycle_Gene[,1])

intra_variance_sum <- a$tot.withinss

wrapk = function(ds, k){
  km = kmeans(ds, k)
  IVsum = km$tot.withinss
  clustered = km$cluster
  clusteredDS = cbind(ds, clustered)
  return(list(clusteredDS, IVsum))
}

x = wrapk(S_cycle, 5)

a = x[[1]]
b = x[[2]]

checkIVsumOf = function(ds, kvec){
  IV = c(1:length(kvec))
  for (i in 1:length(kvec)) {
    clustering = wrapk(ds, kvec[i])
    IV[i] = clustering[[2]]}
  result = cbind(kvec, IV)
  return(result)
}

qwe = checkIVsumOf(S_cycle, c(4,6,8,10,12,14,16))
plot(qwe)

km = kmeans(S_cycle, 4)
clustered = km[[1]]
Heatmap(as.matrix(S_cycle), name = "Heatmap", cluster_rows = FALSE, cluster_columns = FALSE, split = km$cluster,
        show_heatmap_legend = F, gap = unit(5, "mm"))

km_2 = wrapk(S_cycle, 4)
clustered_2 = km_2[[1]]

calculateDistance = function(Data, Method){
  if(Method == "euclidean"){
    outcome = dist(Data, method = Method, diag = FALSE, upper = FALSE, p = 2)
  } else if(Method == "pearson"){
    outcome = as.dist(1 - cor(t(Data), method = Method))
    return(outcome)
  } else if(Method == "spearman"){
    outcome = as.dist(1 - cor(t(Data), method = Method))
    return(outcome)
  } else {
    print("Please enter a correct method name!")
    return(outcome)
  }
}

makeSilPlot = function(ds){
  k = max(ds$cluster)
  original = subset(ds, select = -c(ds$cluster))
  distMat = calculateDistance(original, "euclidean") #��Euclidean������д��������ʽ
  clustering = ds$cluster
  plot(silhouette(clustering, distMat), col = 1:k, border = NA)
}

makeSilPlot(clustered_2)

plot(c(1, 2, 3, 4), 
     c(sum(nrow(subset(clustered_2, clustered_2$clustered == 1))),
                      sum(nrow(subset(clustered_2, clustered_2$clustered == 2))), 
                      sum(nrow(subset(clustered_2, clustered_2$clustered == 3))), 
                      sum(nrow(subset(clustered_2, clustered_2$clustered == 4)))), 
     xlab = "cluster", 
     ylab = "number of genes",
     main = "k = 4, Profile Plot")

# Manufactoring Silhouette plots
makeSilPlot(wrapk(S_cycle, 6)[[1]])
makeSilPlot(wrapk(S_cycle, 8)[[1]])
makeSilPlot(wrapk(S_cycle, 10)[[1]])
makeSilPlot(wrapk(S_cycle, 12)[[1]])
makeSilPlot(wrapk(S_cycle, 14)[[1]])
makeSilPlot(wrapk(S_cycle, 16)[[1]])







