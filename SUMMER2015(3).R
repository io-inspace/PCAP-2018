#SUMMER 2015 
#(3) the HETEROGENEITY-HETEROGENEITY RELATIONSHIP (HHR) hypothesis

#call data etc:
#######################################################################
library(MuMIn)
library(nlme)
library(lme4)
library(RColorBrewer)
library(AICcmodavg)

setwd("~/Documents/Research/PCAP/PCAPbckup7142017/Dataframes")

setwd("~/Documents/Research/PCAP/PCAPdata/Dataframes")
summplot.patterns.a<-read.csv("summplot.patterns.allforms.csv",
                              header=TRUE,sep=",")
summplot.patterns.a<-summplot.patterns.a[,-1]
summplot.patterns.f<-read.csv("summplot.patterns.forbs.csv",header=TRUE,
                              sep=",")
summplot.patterns.f<-summplot.patterns.f[,-1]

summplot.ab.t.z<-read.csv("summplot.ab.t.z.csv",header=TRUE,sep=',')
summplot.ab.t.z<-summplot.ab.t.z[,-1]

#call categorical data reservation and community
setwd("~/Documents/Research/PCAP/PCAPdata")
cats<-read.csv("2015 categorical data.csv",header=TRUE,sep=',')
categories<-matrix(as.numeric(as.factor(as.matrix(cats[,3:4]))),ncol=2)
colnames(categories)<-c("comm","res")

data<-summplot.ab.t.z[,c(18:32,34)]
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:16, 1:16, vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

#a function to run the three diagnostic tests on each model 
diagnostics<- function(x,print=TRUE){
  plot(x) #using Cook's distance to check for outliers
  hist(x$residuals) #to check for normalcy
  a<-coef(x)
  b<-summary(x) #to get coefficients and stats
  c<-AICc(x)
  d<-r.squaredGLMM(x)
  output<-list(a,b,c,d)
  return(output)
}

ctrl<-lmeControl(opt='optim')
options(na.action=na.fail)
#######################################################################
#BETA RICHNESS
#(3)
#######################################################################
#you probably need to decouple these from richness.
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
#need to make sure 10 is removed for all summbeta models
concatenate<-concatenate[c(1:9,11:29),]

#rich<-scale(concatenate$summgamma.rich,center=TRUE)
#data<-data.frame(rich,summplot.ab.t.z[c(1:9,11:29),c(18:32,34)])
#is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
#  if(j >= i) return(NA)
#  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
#  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
#}
#vCorrelated <- Vectorize(is.correlated, c("i", "j"))
#smat <- outer(1:17, 1:17, vCorrelated, data = data)
#smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
#nm <- colnames(data)
#dimnames(smat) <- list(nm, nm)
#richness is correlated with vMg and vCa, so just be aware.

#step (A)
global.model<-lm((1-summbeta.rich) ~ vOM + vP + vKppm + vMgppm + 
                   vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
                   vlight + vl + vr, data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

modA<-lm((1-summbeta.rich)~vMg,data=concatenate)
diagnostics(modA)
#df: 1, 26
#AICc: -38.4
#multiple-r2: 33%

#step (B)
modB<-lme((1-summbeta.rich)~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 24
#AICc: -35.4
#marginal-r2: 0%
#conditional-r2: 38%

#step (C1)
global.model<-lme((1-summbeta.rich) ~ vOM + vP + vKppm + vMgppm + 
                    vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
                    vlight + vl + vr,random=list(comm=pdDiag(~1)),method ='ML', data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

modC1<-lme((1-summbeta.rich) ~ vMg,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modC1)
#df: 1, 23
#AICc: -36.1
#marginal-r2: 21%
#conditional-r2: 31%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[1:17]))){
  predictor<-concatenate[,i+17]
  model<-lme((1-summbeta.rich)~1,random=list(comm=pdDiag(~predictor-1)),control=ctrl,data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,18:34]),score)
View(results[order(results[,2]),])

modC2<-lme((1-summbeta.rich)~1,random=list(comm=pdDiag(~vMg-1)),control=ctrl,data=concatenate,method='ML')
diagnostics(modC2)
#df: 24
#AICc: -35.6
#marginal-r2: 0%
#conditional-r2: 39%

#step (D1)
score<-NULL
for(i in 1:length(colnames(concatenate[1:17]))){
  predictor<-concatenate[,i+17]
  model<-lme((1-summbeta.rich)~vMg,random=list(comm=pdDiag(~predictor)),control=ctrl,data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,18:34]),score)
View(results[order(results[,2]),])

modD1<-lme((1-summbeta.rich) ~ vMg,random=list(comm=pdDiag(~vl)),data=concatenate,method='ML')
diagnostics(modD1)
#df: 1, 23
#AICc: -34.5
#marginal-r2: 6%
#conditional-r2: 51%
modD1.exp<-lme((1-summbeta.rich) ~ vMg+vl,random=list(comm=pdDiag(~vl)),data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 2, 22
#AICc: -31.5
#marginal-r2: 8%
#conditional-r2: 52%

#is vl correlated with the mean of something?
score<-NULL
for(i in 1:length(colnames(concatenate[1:17]))){
  predictor<-concatenate[,i]
  model<-lme((1-summbeta.rich)~vMg,random=list(comm=pdDiag(~predictor)),control=ctrl,data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,1:17]),score)
View(results[order(results[,2]),])

modD1_mean<-lme((1-summbeta.rich) ~ vMg,random=list(comm=pdDiag(~mr)),data=concatenate,method='ML')
diagnostics(modD1_mean)
#df: 1, 23
#AICc: -38.9
#marginal-r2: 22%
#conditional-r2: 53%

#step (D2)
#do this iteratively:
modD2<-lme((1-summbeta.rich)~vOM,random=list(comm=pdDiag(~vMg)),data=concatenate,method='ML')

global.model<-lme((1-summbeta.rich) ~ vOM + vP + vKppm + vMgppm + 
                    vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
                    vlight + vl + vr,random=list(comm=pdDiag(~vMg)),control=ctrl,method ='ML',data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

modD2<-lme((1-summbeta.rich~1),random=list(comm=pdDiag(~vMg)),data=concatenate,method='ML')
diagnostics(modD2)
#df: 24
#AICc: -33.7
#marginal-r2: 0%
#conditional-r2: 37%

modD2.exp<-lme((1-summbeta.rich~vMg),random=list(comm=pdDiag(~vMg)),data=concatenate,method='ML')
diagnostics(modD2.exp)
#df: 1, 23
#AICc: -33.2
#marginal-r2: 22%
#conditional-r2: 34%

#######################################################################
#plot best model(s):
#######################################################################
#best linear model:
plot((1-summbeta.rich)~vMg,data=concatenate,col='black',pch=16,cex=1.7,
     main='Summer HHR model A',
     xlab='normalized magnesium heterogeneity',
     ylab='spatial turnover of species',ylim=c(0,1),cex.axis=1.5)
abline(coef(modA)[1],coef(modA)[2],lwd=5,cex=1.5)
legend('topleft',legend=c('across plots'),lwd=3,cex=1.5,pch=16)

#best lme model:
plot((1-summbeta.rich)~vMg,data=concatenate,pch=NA,
     main='Summer Forbs: the HHR hypothesis 
     within forests',xlab='normalized magnesium heterogeneity',
     ylab = 'spatial turnover of species',cex = 1.2,ylim=c(0,1))
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
abline(0.4850591 ,0.0571163,col='red',lwd=3)
#BM
points(concatenate[concatenate$comm==1,26],
       (1-concatenate[concatenate$comm==1,39]),pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modC1)[1,1],coef(modC1)[1,2],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,26],
       (1-concatenate[concatenate$comm==3,39]),pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modC1)[2,1],coef(modC1)[2,2],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,26],
       (1-concatenate[concatenate$comm==5,39]),pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modC1)[3,1],coef(modC1)[3,2],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,26],
       (1-concatenate[concatenate$comm==8,39]),pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modC1)[4,1],coef(modC1)[4,2],col=colors[4],lwd=3,lty=4)
legend("topleft",legend=c('across plots','Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,1,2,3,4),pch=c(NA,as.character('B'),as.character('F'),as.character('M'),
                            as.character('O')),col=c('red',colors),lwd = c(3,1.5,1.5,1.5,1.5))

#second best model
plot((1-summbeta.rich)~vMg,data=concatenate,pch=NA,
     main='Summer Forbs: the HHR hypothesis 
     within forests',xlab='normalized magnesium heterogeneity',
     ylab = 'spatial turnover of species',cex = 1.2,ylim=c(0,1))
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
abline(0.4883740,0.0316268,col='red',lwd=3)
#BM
points(concatenate[concatenate$comm==1,26],
       (1-concatenate[concatenate$comm==1,39]),pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1)[1,1],coef(modD1)[1,2],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,26],
       (1-concatenate[concatenate$comm==3,39]),pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1)[2,1],coef(modD1)[2,2],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,26],
       (1-concatenate[concatenate$comm==5,39]),pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1)[3,1],coef(modD1)[3,2],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,26],
       (1-concatenate[concatenate$comm==8,39]),pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1)[4,1],coef(modD1)[4,2],col=colors[4],lwd=3,lty=4)
legend("topleft",legend=c('across plots','Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,1,2,3,4),pch=c(NA,as.character('B'),as.character('F'),as.character('M'),
                              as.character('O')),col=c('red',colors),lwd = c(3,1.5,1.5,1.5,1.5))

plot((1-summbeta.rich)~vl,data=concatenate,pch=NA,
     main='Summer Forbs: the HHR hypothesis 
     within forests',xlab='normalized litter depth heterogeneity',
     ylab = 'spatial turnover of species',cex = 1.7,ylim=c(0,1),cex.axis=1.5)
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,32],
       (1-concatenate[concatenate$comm==1,39]),pch=as.character('B'),col=colors[1],cex=1.7)
abline(coef(modD1)[1,1],coef(modD1)[1,3],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,32],
       (1-concatenate[concatenate$comm==3,39]),pch=as.character('F'),col=colors[2],cex=1.7)
abline(coef(modD1)[2,1],coef(modD1)[2,3],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,32],
       (1-concatenate[concatenate$comm==5,39]),pch=as.character('M'),col=colors[3],cex=1.7)
abline(coef(modD1)[3,1],coef(modD1)[3,3],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,32],
       (1-concatenate[concatenate$comm==8,39]),pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modD1)[4,1],coef(modD1)[4,3],col=colors[4],lwd=5,lty=4)
legend("topright",legend=c('Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,2,3,4),pch=c(as.character('B'),as.character('F'),as.character('M'),
                            as.character('O')),col=colors,lwd = c(3,3,3,3),cex=1.5)

#modD1_mean
plot((1-summbeta.rich)~vMg,data=concatenate,pch=NA,
     main='Summer Forbs: the HHR hypothesis 
     within forests',xlab='normalized magnesium heterogeneity',
     ylab = 'spatial turnover of species',cex = 1.2,ylim=c(0,1))
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
abline( 0.4824350,0.0594235,col='red',lwd=3)
#BM
points(concatenate[concatenate$comm==1,26],
       (1-concatenate[concatenate$comm==1,39]),pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1_mean)[1,1],coef(modD1_mean)[1,2],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,26],
       (1-concatenate[concatenate$comm==3,39]),pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1_mean)[2,1],coef(modD1_mean)[2,2],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,26],
       (1-concatenate[concatenate$comm==5,39]),pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1_mean)[3,1],coef(modD1_mean)[3,2],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,26],
       (1-concatenate[concatenate$comm==8,39]),pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1_mean)[4,1],coef(modD1_mean)[4,2],col=colors[4],lwd=3,lty=4)
legend("topleft",legend=c('across plots','Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,1,2,3,4),pch=c(NA,as.character('B'),as.character('F'),as.character('M'),
                              as.character('O')),col=c('red',colors),lwd = c(3,1.5,1.5,1.5,1.5))

plot((1-summbeta.rich)~mr,data=concatenate,pch=NA,
     main='Summer Forbs: the HHR hypothesis 
     within forests',xlab='normalized mean restrictive layer depth',
     ylab = 'spatial turnover of species',cex = 1.2,ylim=c(0,1))
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,17],
       (1-concatenate[concatenate$comm==1,39]),pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1_mean)[1,1],coef(modD1_mean)[1,3],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,17],
       (1-concatenate[concatenate$comm==3,39]),pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1_mean)[2,1],coef(modD1_mean)[2,3],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,17],
       (1-concatenate[concatenate$comm==5,39]),pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1_mean)[3,1],coef(modD1_mean)[3,3],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,17],
       (1-concatenate[concatenate$comm==8,39]),pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1_mean)[4,1],coef(modD1_mean)[4,3],col=colors[4],lwd=3,lty=4)
legend("topright",legend=c('Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,2,3,4),pch=c(as.character('B'),as.character('F'),as.character('M'),
                            as.character('O')),col=colors,lwd = c(1.5,1.5,1.5,1.5))

#######################################################################
#RICHNESS HETEROGENEITY
#(3)
#######################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
concatenate<-concatenate[c(1:23,25:29),]

#rich<-scale(concatenate$summgamma.rich,center=TRUE)
#data<-data.frame(rich,summplot.ab.t.z[c(1:23,25:29),c(18:32,34)])
#is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
#  if(j >= i) return(NA)
#  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
#  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
#}
#vCorrelated <- Vectorize(is.correlated, c("i", "j"))
#smat <- outer(1:17, 1:17, vCorrelated, data = data)
#smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
#nm <- colnames(data)
#dimnames(smat) <- list(nm, nm)
#richness is correlated with vMg and vCa, so just be aware.

#step (A)
global.model<-lm(sqrt(summvar.rich) ~ vOM + vP + vKppm + vMgppm + 
                   vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
                   vlight + vl + vr, data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

modA<-lm(sqrt(summvar.rich)~vMg,data=concatenate) 
diagnostics(modA)
#df: 1, 26
#AICc: 126.6
#multiple-r2: 28%

#step (B)
modB<-lme(sqrt(summvar.rich)~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 24
#AICc: 132.3
#marginal-r2: 0%
#conditional-r2: 29%

#step (C1)
global.model<-lme(sqrt(summvar.rich) ~ vOM + vP + vKppm + vMgppm + 
                    vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
                    vlight + vl + vr,random=list(comm=pdDiag(~1)),method ='ML', data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

modC1<-lme(sqrt(summvar.rich)~vOM+vKppm+vMg,random=list(comm=pdDiag(~1)),method ='ML', data = concatenate)
diagnostics(modC1)
#df: 3, 21
#AICc: 129.0
#marginal-r2: 39%
#conditional-r2: 51%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[18:34]))){
  predictor<-concatenate[,i+17]
  model<-lme(sqrt(summvar.rich)~1,random=list(comm=pdDiag(~predictor-1)),control=ctrl,data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,18:34]),score)
View(results[order(results[,2]),])

modC2<-lme(sqrt(summvar.rich)~1,random=list(comm=pdDiag(~vCEC-1)),control=ctrl,data=concatenate,method='ML')
diagnostics(modC2)
#df: 24
#AICc: 132.3
#marginal-r2: 0%
#conditional-r2: 30%

#step (D1)
score<-NULL
for(i in 1:length(colnames(concatenate[18:34]))){
  predictor<-concatenate[,i+17]
  model<-lme(sqrt(summvar.rich)~vMg+vOM+vKppm,random=list(comm=pdDiag(~predictor)),control=ctrl,data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,18:34]),score)
View(results[order(results[,2]),])

modD1<-lme(sqrt(summvar.rich)~vMg+vOM+vKppm,random=list(comm=pdDiag(~vl)),control=ctrl,data=concatenate,method='ML')
diagnostics(modD1)
#df: 3, 21
#AICc: 130.6
#marginal-r2: 34%
#conditional-r2: 60%
modD1.exp<-lme(sqrt(summvar.rich)~vMg+vOM+vKppm+vl,random=list(comm=pdDiag(~vl)),control=ctrl,data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 4, 20
#AICc: 134.5
#marginal-r2: 34%
#conditional-r2: 60%

#step (D2)
global.model<-lme(sqrt(summvar.rich) ~ vOM + vP + vKppm + vMgppm + 
                    vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
                    vlight + vl + vr,random=list(comm=pdDiag(~vCEC)),method ='ML', data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

modD2<-lme(sqrt(summvar.rich)~vKppm,random=list(comm=pdDiag(~vCEC)),method ='ML', data = concatenate)
diagnostics(modD2)
#df: 1, 23
#AICc: 127.3
#marginal-r2: 11%
#conditional-r2: 74%
modD2.exp<-lme(sqrt(summvar.rich)~vKppm+vCEC,random=list(comm=pdDiag(~vCEC)),method ='ML', data = concatenate)
diagnostics(modD2.exp)
#df: 2, 22
#AICc: 129.1
#marginal-r2: 22%
#conditional-r2: 73%
#######################################################################
#plot best model(s):
#######################################################################
#best model is the linear model:
plot(sqrt(summvar.rich)~vMg,data=concatenate,col='black',pch=16,cex=1.7,
     main='HHR model A',xlab='normalized magnesium heterogeneity',
     ylab='richness heterogeneity',axes=FALSE)
axis(1,at=c(-2,-1,0,1,2),labels=c(-2,-1,0,1,2),cex.axis=1.5)
axis(2,at=c(0,2,4,6,8,10),labels=round((c(0,2,4,6,8,10))^2,1),cex.axis=1.5)
box()
abline(coef(modA)[1],coef(modA)[2],lwd=5,col='black')
legend('bottomright',legend=c('across plots'),lwd=3,col='black',pch=16,cex=1.5)

#second best model is modD2
plot(sqrt(summvar.rich)~vKppm,data=concatenate,pch=NA,
     main='HHR model D2',xlab='normalized potassium (ppm) heterogeneity',
     ylab = 'richness heterogeneity',cex = 1.2,axes=FALSE)
axis(1,at=c(-2,-1,0,1,2,3,4,5),labels=c(-2,-1,0,1,2,3,4,5),cex.axis=1.5)
axis(2,at=c(0,2,4,6,8,10),labels=round((c(0,2,4,6,8,10))^2,1),cex.axis=1.5)
box()
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
abline(5.491752, -0.902320,col='red',lwd=5)
#BM
points(concatenate[concatenate$comm==1,20],
       sqrt(concatenate[concatenate$comm==1,40]),pch=as.character('B'),col=colors[1],cex=1.7)

abline(coef(modD2)[1,1],coef(modD2)[1,2],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,20],
       sqrt(concatenate[concatenate$comm==3,40]),pch=as.character('F'),col=colors[2],cex=1.7)

abline(coef(modD2)[2,1],coef(modD2)[2,2],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,20],
       sqrt(concatenate[concatenate$comm==5,40]),pch=as.character('M'),col=colors[3],cex=1.7)

abline(coef(modD2)[3,1],coef(modD2)[3,2],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,20],
       sqrt(concatenate[concatenate$comm==8,40]),pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modD2)[4,1],coef(modD2)[4,2],col=colors[4],lwd=5,lty=4)
legend("topright",legend=c('across plots','Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,1,2,3,4),pch=c(NA,as.character('B'),as.character('F'),as.character('M'),
                              as.character('O')),col=c('red',colors),lwd = c(3,3,3,3,3),cex=1.5,
       inset=c(-.0175,-.04))

plot(sqrt(summvar.rich)~vCEC,data=concatenate,pch=NA,
     main='HHR model D2',xlab='normalized cation exchange capacity heterogeneity',
     ylab = 'richness heterogeneity',axes=FALSE)
axis(1,at=c(-1.5,-1,-.5,0,.5,1,1.5,2),labels=c(-1.5,-1,-.5,0,.5,1,1.5,2),cex.axis=1.5)
axis(2,at=c(0,2,4,6,8,10),labels=round((c(0,2,4,6,8,10))^2,1),cex.axis=1.5)
box()
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,24],
       sqrt(concatenate[concatenate$comm==1,40]),pch=as.character('B'),col=colors[1],cex=1.7)

abline(coef(modD2)[1,1],coef(modD2)[1,3],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,24],
       sqrt(concatenate[concatenate$comm==3,40]),pch=as.character('F'),col=colors[2],cex=1.7)

abline(coef(modD2)[2,1],coef(modD2)[2,3],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,24],
       sqrt(concatenate[concatenate$comm==5,40]),pch=as.character('M'),col=colors[3],cex=1.7)

abline(coef(modD2)[3,1],coef(modD2)[3,3],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,24],
       sqrt(concatenate[concatenate$comm==8,40]),pch=as.character('O'),col=colors[4],cex=1.7)

abline(coef(modD2)[4,1],coef(modD2)[4,3],col=colors[4],lwd=3,lty=4)
legend("bottomright",legend=c('Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,2,3,4),pch=c(as.character('B'),as.character('F'),as.character('M'),
                              as.character('O')),col=c(colors),lwd = c(1.5,1.5,1.5,1.5))

