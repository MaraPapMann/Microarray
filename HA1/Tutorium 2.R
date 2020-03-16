#Aufgabe a

oligo1 = read.table("C:/Users/isoch/Desktop/OligoHitsums1.csv",header=FALSE,sep=",",row.names=1)

#For Notes
dim(oligo1) #求行，列数量
oligo1[1,] #取第一行
oligo1[,1] #取第一列
oligo1[1,1] #取第一行第一列
oligo1[,1:10] #取所有行，一到十列
oligo1[70:100,] #取第七十行到一百行，所有列
oligo1[,10:30] #取所有行，十到三十列
oligo1[70:100,10:30] #取七十到一百行，十到三十列

#Plotting for the whole
boxplot(oligo1,
        col='green',
        pars=list(outcol='red'),
        main='All probes',
        xlab='number of mismatches',
        ylab='number of oligos with x mismatches')



#Aufgabe b

#Plotting for the probes with 60-100 mismatches
boxplot(oligo1[,59:101],
        col='green',
        pars=list(outcol='red'),
        main='All Probes, 58 to 100 mismatches',
        xlab='number of mismatches',
        ylab='number of oligos with x mismatches')

#As a matter of fact, the number of oligos with x mismatches (0<x<60) is around 0. 
#But the number of oligos with x mismatches (60<= x <=90) are significantly higher.
#When x=75 , the number reaches the highest value at around 1e+05.
#Therefore, for this microarray, a meaningful threshold of mismatches shall be [60,90], which contains the most oligos.


