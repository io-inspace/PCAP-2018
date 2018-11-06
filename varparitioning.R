spring<-c(39,8,28,25)
summer<-c(39,0,51,11)
varpartition<-as.matrix(cbind(spring,summer))
barplot(varpartitioning,col=c('black','slategrey','lightgrey','white'),legend=c('unexplained','heterogeneity','availability','both'))

#SPRING
#across plots: mCa, vP, vK
chart.Correlation(cbind(concatenate[,c(10)],concatenate[,c(19,25)]))
plot(varpart(concatenate$spgamma.rich,concatenate[,c(10)],concatenate[,c(19,25)],data=concatenate))
#28, 25, 8

#within forests: 
#Beech-Maple
concatenateBM<-concatenate[concatenate$comm==1,]
chart.Correlation(cbind(concatenateBM[,c(10,12,17)],concatenateBM[,c(19,25)]))
#MINIMIZING RESIDUALS
plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(19)],data=concatenateBM))
plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(25)],data=concatenateBM))
plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(12)],concatenateBM[,c(19)],data=concatenateBM))
plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(12)],concatenateBM[,c(25)],data=concatenateBM))
plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(17)],concatenateBM[,c(19)],data=concatenateBM))
plot(varpart(concatenateBM$spgamma.rich,concatenateBM[,c(17)],concatenateBM[,c(25)],data=concatenateBM))

concatenateF<-concatenate[concatenate$comm==3,]
chart.Correlation(cbind(concatenateF[,c(10,12,17)],concatenateF[,c(19,25)]))

plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(10)],concatenateF[,c(19)],data=concatenateF))
#53
plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(10)],concatenateF[,c(25)],data=concatenateF))
#38
plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(12)],concatenateF[,c(19)],data=concatenateF))
#61
plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(12)],concatenateF[,c(25)],data=concatenateF))
#46
plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(17)],concatenateF[,c(19)],data=concatenateF))
#52
plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(17)],concatenateF[,c(25)],data=concatenateF))
#35

#plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(10)],concatenateF[,c(19,25)],data=concatenateF))
##3, 46, 8
#plot(varpart(concatenateF$spgamma.rich,concatenateF$mCa,concatenateF$vP,concatenateF$vK,data=concatenateF))
#10 - vK
#plot(varpart(concatenateF$spgamma.rich,concatenateF$mN,concatenateF$mr,data=concatenateF))
#0, 0, 40
#plot(varpart(concatenateF$spgamma.rich,concatenateF[,c(17)],concatenateF[,c(19,25)],data=concatenateF))
#6, 31, 23

concatenateM<-concatenate[concatenate$comm==5,]
chart.Correlation(cbind(concatenateM[,c(10,12,17)],concatenateM[,c(19,25)]))

plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(10)],concatenateM[,c(19)],data=concatenateM))
#1.67
plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(10)],concatenateM[,c(25)],data=concatenateM))
#1.44
plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(12)],concatenateM[,c(19)],data=concatenateM))
#.15
plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(12)],concatenateM[,c(25)],data=concatenateM))
#.01, m>1
plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(17)],concatenateM[,c(19)],data=concatenateM))
#.21
plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(17)],concatenateM[,c(25)],data=concatenateM))
#0, m>1

concatenateO<-concatenate[concatenate$comm==8,]
chart.Correlation(cbind(concatenateO[,c(10,12,17)],concatenateO[,c(19,25)]))

plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(10)],concatenateO[,c(19)],data=concatenateO))
#.4
plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(10)],concatenateO[,c(25)],data=concatenateO))
#.63
plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(12)],concatenateO[,c(19)],data=concatenateO))
#.50
plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(12)],concatenateO[,c(25)],data=concatenateO))
#.54
plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(17)],concatenateO[,c(19)],data=concatenateO))
#.53
plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(17)],concatenateO[,c(25)],data=concatenateO))
#.58

sp.bm<-c(18,80,29,0)
sp.fp<-c(35,29,5,32)
sp.m<-c(15,0,83,10)
sp.oak<-c(40,15,27,18)
sp.varpartition<-as.matrix(cbind(sp.bm,sp.fp,sp.m,sp.oak))
barplot(sp.varpartition,col=c('black','slategrey','lightgrey','white'))#,legend=c('unexplained','heterogeneity','availability','both'))
 

#plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(10)],concatenateM[,c(19,25)],data=concatenateM))
#plot(varpart(concatenateM$spgamma.rich,concatenateM$mCa,concatenateM$vP,concatenateM$vK,data=concatenateM))
#plot(varpart(concatenateM$spgamma.rich,concatenateM$mN,concatenateM$mr,data=concatenateM))
#17, 76, 7
#plot(varpart(concatenateM$spgamma.rich,concatenateM[,c(12)],concatenateM[,c(19)],data=concatenateM))

#concatenateO<-concatenate[concatenate$comm==8,]
#plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(10)],concatenateO[,c(19,25)],data=concatenateO))

#plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(10,12,17)],concatenateO[,c(19,25)],data=concatenateO))
#29, 16, 7
#plot(varpart(concatenateO$spgamma.rich,concatenateO$mCa,concatenateO$vP,concatenateO$vK,data=concatenateO))
#29 - mCa
#15 - vP
#plot(varpart(concatenateO$spgamma.rich,concatenateO$mN,concatenateO$mr,data=concatenateO))
#50, 4, 23
#plot(varpart(concatenateO$spgamma.rich,concatenateO[,c(12,17)],concatenateO[,c(19,25)],data=concatenateO))
#66,11,12
#fit = lmer(spgamma.rich~ mCa+vP+(mr | comm),data=concatenate)
#VarCorrCI(fit)
#???

#SUMMER
#plot-level
chart.Correlation(cbind(concatenate[,c(2,8,10)],concatenate[,c(23,29)]))
plot(varpart(concatenate$summgamma.rich,concatenate[,c(2,8,10)],concatenate[,c(23,29)],data=concatenate))

concatenateBM<-concatenate[concatenate$comm==1,]
chart.Correlation(cbind(concatenateBM[,c(2,8,10)],concatenateBM[,c(23,24,29)]))

plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2)],concatenateBM[,c(23)],data=concatenateBM))
#.59
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2)],concatenateBM[,c(24)],data=concatenateBM))
#1.44
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2)],concatenateBM[,c(29)],data=concatenateBM))
#.96
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(8)],concatenateBM[,c(23)],data=concatenateBM))
#.65
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(8)],concatenateBM[,c(24)],data=concatenateBM))
#1.22
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(8)],concatenateBM[,c(29)],data=concatenateBM))
#.93
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(23)],data=concatenateBM))
#.84
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(24)],data=concatenateBM))
#1.28
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(29)],data=concatenateBM))
#1.05

concatenateF<-concatenate[concatenate$comm==3,]
chart.Correlation(cbind(concatenateF[,c(2,8,10)],concatenateF[,c(23,24,29)]))

plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(2)],concatenateF[,c(23)],data=concatenateF))
#.81
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(2)],concatenateF[,c(24)],data=concatenateF))
#.69
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(2)],concatenateF[,c(29)],data=concatenateF))
#.65
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(8)],concatenateF[,c(23)],data=concatenateF))
#.40
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(8)],concatenateF[,c(24)],data=concatenateF))
#.67
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(8)],concatenateF[,c(29)],data=concatenateF))
#.70
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(10)],concatenateF[,c(23)],data=concatenateF))
#.72
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(10)],concatenateF[,c(24)],data=concatenateF))
#1.02
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(10)],concatenateF[,c(29)],data=concatenateF))
#.86


concatenateM<-concatenate[concatenate$comm==5,]
chart.Correlation(cbind(concatenateM[,c(2,8,10)],concatenateM[,c(23,24,29)]))

plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2)],concatenateM[,c(23)],data=concatenateM))
#1.51
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2)],concatenateM[,c(24)],data=concatenateM))
#.86
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2)],concatenateM[,c(29)],data=concatenateM))
#.56
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(8)],concatenateM[,c(23)],data=concatenateM))
#1.34
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(8)],concatenateM[,c(24)],data=concatenateM))
#.87
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(8)],concatenateM[,c(29)],data=concatenateM))
#.59
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(10)],concatenateM[,c(23)],data=concatenateM))
#1.43
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(10)],concatenateM[,c(24)],data=concatenateM))
#.36
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(10)],concatenateM[,c(29)],data=concatenateM))
#.57

concatenateO<-concatenate[concatenate$comm==8,]
chart.Correlation(cbind(concatenateO[,c(2,8,10)],concatenateO[,c(23,24,29)]))

plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2)],concatenateO[,c(23)],data=concatenateO))
#.32
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2)],concatenateO[,c(24)],data=concatenateO))
#1.29
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2)],concatenateO[,c(29)],data=concatenateO))
#1.08
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(8)],concatenateO[,c(23)],data=concatenateO))
#.30
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(8)],concatenateO[,c(24)],data=concatenateO))
#1.30
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(8)],concatenateO[,c(29)],data=concatenateO))
#1.08
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(10)],concatenateO[,c(23)],data=concatenateO))
#.32
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(10)],concatenateO[,c(24)],data=concatenateO))
#.69
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(10)],concatenateO[,c(29)],data=concatenateO))
#.63

summ.bm<-c(59,49,4,0)
summ.fp<-c(40,22,31,7)
summ.m<-c(36,96,26,0)
summ.oak<-c(30,82,0,0)
summ.varpartition<-as.matrix(cbind(summ.bm,summ.fp,summ.m,summ.oak))
barplot(summ.varpartition,col=c('black','slategrey','lightgrey','white'))#,legend=c('unexplained','heterogeneity','availability','both'))





#51,11,0
plot(varpart(concatenate$summgamma.rich,concatenate[,c(8,10)],concatenate[,c(23,29)],data=concatenate))
#49, 5, 5
plot(varpart(concatenate$summgamma.rich,concatenate[,c(2,8)],concatenate[,c(23,29)],data=concatenate))
#3, 0, 16
plot(varpart(concatenate$summgamma.rich,concatenate$mCa,concatenate$mK,concatenate$vpH,concatenate$vN,data=concatenate))
plot(varpart(concatenate$summgamma.rich,concatenate$mP,concatenate$mK,concatenate$vpH,concatenate$vN,data=concatenate))
#COMPARABLE
plot(varpart(concatenate$summgamma.rich,concatenate[,c(10)],concatenate[,c(23,29)],data=concatenate))
#38,8,2
#plot(varpart(concatenate$summgamma.rich,concatenate$mN,concatenate$mr,data=concatenate))


concatenateBM<-concatenate[concatenate$comm==1,]
#collinearity: mCa and mK, vpH and vN

#mCa and mK 
#vpH and vN 
#are correlated in BM forests
plot(varpart(concatenateBM$summgamma.rich,concatenateBM$mCa,concatenateBM$mP,concatenateBM$vpH,concatenateBM$vN,data=concatenateBM))
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2,8,10)],concatenateBM[,c(23,24,29)],data=concatenateBM))
#choose mK and vpH
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2,8)],concatenateBM[,c(23,24)],data=concatenateBM))
#33, 0, 70, 50
#choose mCa and vpH
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2,10)],concatenateBM[,c(23,24)],data=concatenateBM))
#34, 0, 53, 49
#plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2,8)],concatenateBM[,c(24,29)],data=concatenateBM))
#plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(2,10)],concatenateBM[,c(24,29)],data=concatenateBM))

#COMPARABLE
plot(varpart(concatenateBM$summgamma.rich,concatenateBM[,c(10)],concatenateBM[,c(23,29)],data=concatenateBM))

concatenateF<-concatenate[concatenate$comm==3,]
plot(varpart(concatenateF$summgamma.rich,concatenateF$mCa,concatenateF$mK,concatenateF$vpH,concatenateF$vN,data=concatenateF))
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(8,10)],concatenateF[,c(23,24,29)],data=concatenateF))
#10, 58,0
#COMPARABLE
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(10)],concatenateF[,c(23,29)],data=concatenateF))

concatenateM<-concatenate[concatenate$comm==5,]
#plot(varpart(concatenateM$summgamma.rich,concatenateM$mCa,concatenateM$mP,concatenateM$vpH,concatenateM$vN,data=concatenateM))
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(10)],concatenateM[,c(23,29)],data=concatenateM))
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2)],concatenateM[,c(23,29)],data=concatenateM))
#0, 0, 89
#plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(8)],concatenateM[,c(23,29)],data=concatenateM))

plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(10)],concatenateM[,c(23,29)],data=concatenateM))
#0, 0, 94
#plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2,8)],concatenateM[,c(23,29)],data=concatenateM))

#plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2,10)],concatenateM[,c(23,29)],data=concatenateM))
#plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(8,10)],concatenateM[,c(23,29)],data=concatenateM))

concatenateO<-concatenate[concatenate$comm==8,]
plot(varpart(concatenateO$summgamma.rich,concatenateO$mCa,concatenateO$mP,concatenateO$vpH,concatenateO$vN,data=concatenateO))
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2,8,10)],concatenateO[,c(23,24,29)],data=concatenateO))
#0, 73, 0
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2)],concatenateO[,c(23,29)],data=concatenateO))
#0, 0, 79
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(8)],concatenateO[,c(23,29)],data=concatenateO))
#0, 0, 83
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(10)],concatenateO[,c(23,29)],data=concatenateO))
#0, 45, 29
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2,8)],concatenateO[,c(23,29)],data=concatenateO))
#0, 0, 92
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2,10)],concatenateO[,c(23,29)],data=concatenateO))
#0, 43, 31
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(8,10)],concatenateO[,c(23,29)],data=concatenateO))
#0, 65, 5






concatenateF<-concatenate[concatenate$comm==3,]
plot(varpart(concatenateF$summgamma.rich,concatenateF$mCa,concatenateF$mP,concatenateF$vpH,concatenateF$vN,data=concatenateF))
plot(varpart(concatenateF$summgamma.rich,concatenateF[,c(2,8,10)],concatenateF[,c(23,29)],data=concatenateF))
concatenateM<-concatenate[concatenate$comm==5,]
plot(varpart(concatenateM$summgamma.rich,concatenateM$mCa,concatenateM$mP,concatenateM$vpH,concatenateM$vN,data=concatenateM))
plot(varpart(concatenateM$summgamma.rich,concatenateM[,c(2,8,10)],concatenateM[,c(23,29)],data=concatenateM))
concatenateO<-concatenate[concatenate$comm==8,]
plot(varpart(concatenateO$summgamma.rich,concatenateO$mCa,concatenateO$mP,concatenateO$vpH,concatenateO$vN,data=concatenateO))
plot(varpart(concatenateO$summgamma.rich,concatenateO[,c(2,8,10)],concatenateO[,c(23,29)],data=concatenateO))

fit = lmer(spgamma.rich~ mCa+vP+(mr | comm),data=concatenate)
VarCorrCI(fit)
#???