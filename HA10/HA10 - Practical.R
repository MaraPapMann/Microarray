# FDR 

# Data Reading
setwd("C:/Users/isoch/Google 云端硬盘/Study/1. Semester an Uni-Tuebingen/Microarray Bioinformatik/HA10")
Original_Data <- read.table("golub.dat")
Colname123 <- read.table("colnames.txt")
colnames(Original_Data) <- Colname123[,1]

# Parameters for t-test ()

# 把每一行的ALL和AML视为两个独立不相关的分布，假设他们符合同一个分布，然后进行T检验，不过这里假设他们方差相同。
# 所以进行的是Student t-test.

# t.test(Orginal_Data[1, 1:27], Original_Data[1, 28:38], var.equal = TRUE, "less", conf.level = 0.999) Example

ALL <- Original_Data[,1:27]
AML <- Original_Data[,28:38]

# Programming t-test

self_made_t_test = function(x, y, kind_of_test, alpha){
  base = NULL
  for (i in c(1:nrow(x))) {
    test_outcome = t.test(x[i,], y[i,], var.equal = T, kind_of_test, conf.level = 1-alpha)
    p = test_outcome$p.value
    if (kind_of_test == "two.sided"){
      if (p<(alpha/2) | p>(1-alpha/2)) {
        rejection = T
      } else {rejection = F}
    } else if (kind_of_test == "greater"){
      if (p<alpha){
        rejection = T
      } else {rejection = F}
    } else if (kind_of_test == "less"){
      if (p<alpha){
        rejection = T
      } else {rejection = F}
    }
    base = rbind(base, rejection)
  }
  return(base)
}

two_sided_test_1 <- self_made_t_test(ALL, AML, "two.sided", 0.05)
two_sided_test_2 <- self_made_t_test(ALL, AML, "two.sided", 0.01)
two_sided_test_3 <- self_made_t_test(ALL, AML, "two.sided", 0.001)

Number_of_differentially_expressed_genes1 <- sum(two_sided_test_1 == T)
Number_of_differentially_expressed_genes2 <- sum(two_sided_test_2 == T)
Number_of_differentially_expressed_genes3 <- sum(two_sided_test_3 == T)

greater_test_1 <- self_made_t_test(ALL, AML, "greater", 0.05)
greater_test_2 <- self_made_t_test(ALL, AML, "greater", 0.01)
greater_test_3 <- self_made_t_test(ALL, AML, "greater", 0.001)

Number_of_up_regulated_genes1 <- sum(greater_test_1 == T)
Number_of_up_regulated_genes2 <- sum(greater_test_2 == T)
Number_of_up_regulated_genes3 <- sum(greater_test_3 == T)

less_test_1 <- self_made_t_test(ALL, AML, "less", 0.05)
less_test_2 <- self_made_t_test(ALL, AML, "less", 0.01)
less_test_3 <- self_made_t_test(ALL, AML, "less", 0.001)

Number_of_down_regulated_genes1 <- sum(less_test_1 == T)
Number_of_down_regulated_genes2 <- sum(less_test_2 == T)
Number_of_down_regulated_genes3 <- sum(less_test_3 == T)

# making the table
final_table <- as.matrix(rbind(c(Number_of_differentially_expressed_genes1, 
                                   Number_of_differentially_expressed_genes2, 
                                   Number_of_differentially_expressed_genes3),
                                 c(Number_of_up_regulated_genes1,
                                   Number_of_up_regulated_genes2,
                                   Number_of_up_regulated_genes3),
                                 c(Number_of_down_regulated_genes1,
                                   Number_of_down_regulated_genes2,
                                   Number_of_down_regulated_genes3)))
colnames(final_table) <- c("alpha 0.05", "alpha 0.01", "alpha 0.001")
rownames(final_table) <- c("differentially_expressed_genes", "up_regulated_genes", "down_regulated_genes")

# If alpha = 0.05, the up-regulated genes are more than the down-regulated in ALL. In other cases it is on the contrary.
# The number changes significantly with different alpha.
