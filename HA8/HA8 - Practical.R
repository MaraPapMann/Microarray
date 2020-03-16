# Data Reading
library(gplots)

setwd("C:/Users/isoch/Google ‘∆∂À”≤≈Ã/Study/1. Semester an Uni-Tuebingen/Microarray Bioinformatik/HA8/material")

Spell_Man_Cell_Cycle_Fullset <- read.table("Spellman_cellcycle_alpha_completeexpressionmatrix.txt", header = T)

Name_list_Cell_Cycle_Gene <- read.table("Spellman_cellcycle_alpha_subsetcellcycle_geneIDs.txt")

S_cycle <- subset(Spell_Man_Cell_Cycle_Fullset, row.names(Spell_Man_Cell_Cycle_Fullset) %in% Name_list_Cell_Cycle_Gene[,1])

S_random <- Spell_Man_Cell_Cycle_Fullset[500:(500+nrow(S_cycle)-1),]

# Distance Function Programming

Distance = function(Data, Method){
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

# Calculating

S_c_euclidean_dis <- as.matrix(Distance(S_cycle, "euclidean"))
S_r_euclidean_dis <- as.matrix(Distance(S_random, "euclidean"))

S_c_pearson_dis <- as.matrix(Distance(S_cycle, "pearson"))
S_r_pearson_dis <- as.matrix(Distance(S_random, "pearson"))

S_c_spearman_dis <- as.matrix(Distance(S_cycle, "spearman"))
S_r_spearman_dis <- as.matrix(Distance(S_random, "spearman"))

# QQplot Plotting

qqplot(S_c_euclidean_dis, S_r_euclidean_dis) #QQplot_Euclidean
qqplot(S_c_pearson_dis, S_r_pearson_dis) #QQplot_Pearson
qqplot(S_c_spearman_dis, S_r_spearman_dis) #QQplot_Spearman

# Hierarchical_Clustering Function Programming

Hierarchical_Clustering = function(Data, Cluster_Method, Distance_Method){
  di = Distance(Data, Distance_Method)
  clustered = hclust(di, method = Cluster_Method, members = NULL)
  return(clustered)
}

# Clustering

E_dist_clustered_genome <- Hierarchical_Clustering(S_cycle, "average", "euclidean")
P_dist_clustered_genome <- Hierarchical_Clustering(S_cycle, "average", "pearson")

# Editing Color
library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# Generating Heat Map

# Average/Euclidean/Heatmap
heatmap.2(as.matrix(S_cycle[as.array(E_dist_clustered_genome[["order"]]),]),
          main="Average/Euclidean/Heatmap",
          trace="none",
          col=hmcol)

# Average/Pearson/Heatmap
heatmap.2(as.matrix(S_cycle[as.array(P_dist_clustered_genome[["order"]]),]),
          main="Average/Pearson/Heatmap",
          trace="none",
          col=hmcol)


