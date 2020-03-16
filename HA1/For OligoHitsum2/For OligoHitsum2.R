#Aufgabe a

oligo2 = read.table("C:/Users/isoch/Desktop/OligoHitsums2.csv",header=FALSE,sep=",",row.names=1)
dim(oligo2)

#Plotting for the whole
boxplot(oligo2,
        col='green',
        pars=list(outcol='red'),
        main='All probes',
        xlab='number of mismatches',
        ylab='number of oligos with x mismatches')



#Aufgabe b

#Plotting for the probes with 60-100 mismatches
boxplot(oligo2[,34:61],
        col='green',
        pars=list(outcol='red'),
        main='All Probes, 33 to 58 mismatches',
        xlab='number of mismatches',
        ylab='number of oligos with x mismatches')

#As a matter of fact, the number of oligos with x mismatches (0<x<33) is around 0. 
#But the number of oligos with x mismatches (33<= x <=58) are significantly higher.
#When x=47 , the number reaches the highest value at around 1e+05.
#Therefore, for this microarray, a meaningful threshold of mismatches shall be (33<= x <=58), which contains most of the oligos.
#From the biological point of view, this threshold means that for the probes with the length of 60 in this microarray, the oligos with x mismatches concerntrate within (33,58) and reach its peak at 47.



#Aufgabe c

length(which(oligo2 [,1] == 1))
length(which(oligo2 [,1] > 1))
print('The number of oligos with only one perfect match is 812289')
print('The number of oligos with more than one perfect match is 9451')


barplot(table(oligo2$V2),
        main = 'Oligos with Different Numbers of Perfect Match',
        xlab = 'Numbers of Perfect Match',
        ylab = 'Numbers of Probes',
)

#Aufgabe d
#Sp = |TN|/(|TN|+|FP|) , Sensitivity should also be considered
