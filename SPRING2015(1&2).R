#SPRING 2015 
#combination models (AE & HDR hypotheses)

#call data etc:
########################################################################
#libraries+
library(MuMIn)
library(nlme)
library(lme4)
library(RColorBrewer)
library(AICcmodavg)
library(car)

#call data
setwd("~/Documents/Research/PCAP/PCAPdata/Dataframes")

spplot.patterns.a<-read.csv("spplot.patterns.allforms.csv",
                            header=TRUE,sep=",")
spplot.patterns.a<-spplot.patterns.a[,-1]
spplot.patterns.f<-read.csv("spplot.patterns.forbs.csv",header=TRUE,
                            sep=",")
spplot.patterns.f<-spplot.patterns.f[,-1]
spplot.patterns.o<-read.csv("spplot.patterns.otherforms.csv",header=TRUE,sep=',')
spplot.patterns.o<-spplot.patterns.o[,-1]

spplot.ab.t.z<-read.csv("spplot.ab.t.z.csv",header=TRUE,sep=',')
spplot.ab.t.z<-spplot.ab.t.z[,-1]

#call categorical data reservation and community
setwd("~/Documents/Research/PCAP/PCAPdata")
cats<-read.csv("2015 categorical data.csv",header=TRUE,sep=',')
categories<-matrix(as.numeric(as.factor(as.matrix(cats[,3:4]))),ncol=2)
colnames(categories)<-c("comm","res")

#set up correlation matrix for abiotic factors
data<-spplot.ab.t.z[,c(10,25,19)] #includes mCa, vK, vP
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

#a function to run diagnostic tests on each model 
diagnostics<- function(x,print=TRUE){
  plot(x) #using Cook's distance to check for outliers
  hist(x$residuals) #to check for normalcy
  a<-coef(x)#to get model coefficients
  b<-summary(x)#to get coefficients and stats
  c<-AICc(x)#self-explanatory
  d<-r.squaredGLMM(x)#r2 values for lme models only!  Use r2 stats from summary for lm models. 
  output<-list(a,b,c,d)
  return(output)
}

ctrl<-lmeControl(opt='optim')#ctrl not implemented, just in case.
options(na.action=na.fail)
########################################################################
#GAMMA RICHNESS
#combination model
########################################################################
concatenate<-data.frame(spplot.ab.t.z,categories,spplot.patterns.f)
#remove 10 to match other models so AICc are comparable
concatenate<-concatenate[c(1:9,11:29),]

#step (A)
global.model<-lm(spgamma.rich~mCa+vK+vP,data=concatenate)
model1<-dredge(global.model,subset=smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

modA<-lm(spgamma.rich~mCa+vP,data=concatenate)
diagnostics(modA)
#df: 2, 25
#AICc: 182.1
#adjusted-r2: 56%

#step (B)
#same as in relationships (1) & (2)

#step (C1)
global.model<-lme(spgamma.rich~mCa+vK+vP,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
model1<-dredge(global.model,subset=smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

#the best model is the HDR model, but we're interested in the mean and 
#heterogeneity combined so we report the best combination model, which is:
modC1<-lme(spgamma.rich~mCa+vP,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modC1)
#df: 2, 22
#AICc: 185.1
#marginal-r2: 60%
#conditional-r2: 60%

#step (C2)
#same as in relationships (1) & (2) and the AE model has a lower AICc
score<-NULL
for(i in 1:length(concatenate[,1:34])){
  predictor<-concatenate[,i]
  model<-lme(spgamma.rich~1,random=list(comm=pdDiag(~predictor-1)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,1:34]),score)
View(results[order(results[,2]),])

#best slopes model:
modC2<-lme(spgamma.rich~1,random=list(comm=pdDiag(~mCa-1)),data=concatenate,method='ML')
diagnostics(modC2)
#same as for relationship (1)

#step (D1)
score<-NULL
for(i in 1:length(concatenate[,1:17])){
  predictor<-concatenate[,i]
  #use the best model here, and since it is the HDR model, only let slopes
  #vary by the mean of abiotic factors to get a combination model:
  model<-lme(spgamma.rich~vK+vP,random=list(comm=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,1:17]),score)
View(results[order(results[,2]),])

#best model constrained by intercepts:
modD1<-lme(spgamma.rich~vK+vP,random=list(comm=pdDiag(~mr)),data=concatenate,method='ML')
diagnostics(modD1)
#df: 2, 22
#AICc: 178.2
#marginal-r2: 25%
#conditional-r2: 83%
modD1.exp<-lme(spgamma.rich~vK+vP+mr,random=list(comm=pdDiag(~mr)),data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 3, 21
#AICc: 181.7
#marginal-r2: 27%
#conditional-r2: 83%

#for completeness:
score<-NULL
for(i in 1:length(colnames(concatenate[1:34]))){
  predictor<-concatenate[,i]
  model<-lme(spgamma.rich~mCa+vP,random=list(comm=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}#this formulation of the random parts of the model helps control for correlation between intercepts and slopes.
results<-data.frame(colnames(concatenate[,1:34]),score)
View(results[order(results[,2]),])

modD1.l<-lme(spgamma.rich~mCa+vP,random=list(comm=pdDiag(~vl)),data=concatenate,method='ML')
diagnostics(modD1.l)
#df: 2, 22
#AICc: 184.7
#marginal-r2: 54%
#conditional-r2: 70%
modD1.l.exp<-lme(spgamma.rich~mCa+vP+vl,random=list(comm=pdDiag(~vl)),data=concatenate,method='ML')
diagnostics(modD1.l.exp)
#df: 3, 21
#AICc: 188.1
#marginal-r2: 54%
#conditional-r2: 69%

#step (D2)
global.model<-lme(spgamma.rich~mCa+vP+vK,random=list(comm=pdDiag(~mCa)),data=concatenate,method='ML')
model1<-dredge(global.model,subset=smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

#best model constrained by slopes:
modD2<-lme(spgamma.rich~vK+vP,random=list(comm=pdDiag(~mCa)),data=concatenate,method='ML')
diagnostics(modD2)
#df: 2, 22
#AICc: 186.5
#marginal-r2: 25%
#conditional-r2: 69%
modD2.exp<-lme(spgamma.rich~mCa+vK+vP,random=list(comm=pdDiag(~mCa)),data=concatenate,method='ML')
diagnostics(modD2.exp)
#df: 3, 21
#AICc: 185.9
#marginal-r2: 50%
#conditional-r2: 68%
########################################################################
#plot best model(s):
########################################################################
#modA
plot(spgamma.rich~mCa,data=concatenate,pch=16,cex=1.7,
     main='combination model (A)',xlab='normalized mean calcium',
     ylab = 'plot-level richness (S)',ylim=c(0,35),cex.axis=1.5)
abline(coef(modA)[1],coef(modA)[2],lwd=5)
legend('topleft',legend=c('across plots'),lwd=3,cex=1.5,pch=16)

plot(spgamma.rich~vP,data=concatenate,pch=16,cex=1.7,
     main='combination model (A)',xlab='normalized phosphorus heterogeneity',
     ylab = 'plot-level richness (S)')
abline(coef(modA)[1],coef(modA)[3],lwd=5)

#best combination model modD1:
vif(modD1) #all<2
plot(spgamma.rich~vK,data=concatenate,pch=NA,cex=1.2,
     main='combination model (D1)',
     xlab='normalized potassium heterogeneity',
     ylab = 'plot-level richness (S)',xlim=c(-2.5,2.5),ylim=c(0,35),cex.axis=1.5)
abline(14.575060,-3.712780,col='red',lwd=5)
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,25],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.7)
abline(coef(modD1)[1,1],coef(modD1)[1,2],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,25],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.7)
abline(coef(modD1)[2,1],coef(modD1)[2,2],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,25],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.7)
abline(coef(modD1)[3,1],coef(modD1)[3,2],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,25],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modD1)[4,1],coef(modD1)[4,2],col=colors[4],lwd=5,lty=4)
legend("topright",legend=c('across plots','Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,1,2,3,4),pch=c(NA,as.character('B'),as.character('F'),as.character('M'),
                              as.character('O')),cex=1.5,col=c('red',colors),lwd=c(3,3,3,3,3))

plot(spgamma.rich~vP,data=concatenate,pch=NA,cex=1.7,
     main='combination model (D1)',
     xlab='normalized phosphorus heterogeneity',
     ylab = 'plot-level richness (S)',ylim=c(0,35),cex.axis=1.5)
abline(14.575060,-1.144592,col='red',lwd=5)
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,19],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.7)
abline(coef(modD1)[1,1],coef(modD1)[1,3],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,19],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.7)
abline(coef(modD1)[2,1],coef(modD1)[2,3],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,19],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.7)
abline(coef(modD1)[3,1],coef(modD1)[3,3],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,19],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modD1)[4,1],coef(modD1)[4,3],col=colors[4],lwd=5,lty=4)

plot(spgamma.rich~mr,data=concatenate,pch=NA,cex=1.7,
     main='combination model (D1)',xlab='normalized mean restrictive layer depth',
     ylab = 'plot-level richness (S)',ylim=c(0,35),cex.axis=1.5)
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,17],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.7)
abline(coef(modD1)[1,1],coef(modD1)[1,4],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,17],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.7)
abline(coef(modD1)[2,1],coef(modD1)[2,4],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,17],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.7)
abline(coef(modD1)[3,1],coef(modD1)[3,4],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,17],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modD1)[4,1],coef(modD1)[4,4],col=colors[4],lwd=5,lty=4)

########################################################################
#GAMMA EVENNESS
#combination model
########################################################################
concatenate<-data.frame(spplot.ab.t.z,categories,spplot.patterns.f)
#removed 10 because it was an outlier in the HDR model, 24 for NAs
concatenate<-concatenate[c(1:9,11:23,25:29),] 
rich<-scale(concatenate$spgamma.rich,center=TRUE)
#data<-data.frame(rich,spplot.ab.t.z[c(1:9,11:23,25:29),c(2,7,19,32)])

data<-data.frame(spplot.ab.t.z[c(1:9,11:23,25:29),c(2,7,19,32)])
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

#step (A)
global.model<-lm(log(spgamma.evens)~mP+mCEC+vP+vl,data=concatenate)
model1<-dredge(global.model,subset=smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

#this is not the best model from dredge (that would be the HDR model),
#but it's the best linear model including both abiotic quality and heterogeneity 
modA<-lm(log(spgamma.evens)~vP+vl+mCEC,data=concatenate)
diagnostics(modA)
#df: 3, 23
#AICc: 40.9
#adjusted-r2: 47%
########################################################################
#plot best model(s):
########################################################################
vif(modA) #all <2
plot(log(spgamma.evens)~vP,data=concatenate,pch=16,cex=1.2,
     main='Spring Forbs: combination model and evenness 
     across plots',
     xlab='normalized phosphorus heterogeneity',
     ylab='plot-level evenness (E)',axes=FALSE)
axis(1,at=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
axis(2,at=c(-3,-2,-1.5,-1,-.5,0),labels=round(exp(c(-3,-2,-1.5,-1,-.5,0)),1))
box()
abline(coef(modA)[1],coef(modA)[2],lwd=3)

plot(log(spgamma.evens)~vl,data=concatenate,pch=16,cex=1.2,
     main='Spring Forbs: combination model and evenness 
     across plots',
     xlab='normalized litter depth heterogeneity',
     ylab='plot-level evenness (E)',axes=FALSE)
axis(1,at=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
axis(2,at=c(-3,-2,-1.5,-1,-.5,0),labels=round(exp(c(-3,-2,-1.5,-1,-.5,0)),1))
box()
abline(coef(modA)[1],coef(modA)[3],lwd=3)
legend('topright',legend=c('across plots'),pch=16,lwd=1.5,col='black')

plot(log(spgamma.evens)~mCEC,data=concatenate,pch=16,cex=1.2,
     main='Spring Forbs: combination model and evenness 
     across plots',
     xlab='normalized mean cation exchange capacity',
     ylab='plot-level evenness (E)',axes=FALSE)
axis(1,at=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
axis(2,at=c(-3,-2,-1.5,-1,-.5,0),labels=round(exp(c(-3,-2,-1.5,-1,-.5,0)),1))
box()
abline(coef(modA)[1],coef(modA)[4],lwd=3)

