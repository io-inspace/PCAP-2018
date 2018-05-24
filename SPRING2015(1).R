#SPRING 2015
#(1) PLOT-LEVEL AVAILABLE ENERGY (AE) hypothesis

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

#raw community pattern data:
spplot.patterns.a<-read.csv("spplot.patterns.allforms.csv",header=TRUE,sep=',')
spplot.patterns.a<-spplot.patterns.a[,-1]
spplot.patterns.f<-read.csv("spplot.patterns.forbs.csv",header=TRUE,sep=',')
spplot.patterns.f<-spplot.patterns.f[,-1]
spplot.patterns.o<-read.csv("spplot.patterns.otherforms.csv",header=TRUE,sep=',')
spplot.patterns.o<-spplot.patterns.o[,-1]

#transformed/normalized abiotic data:
spplot.ab.t.z<-read.csv("spplot.ab.t.z.csv",header=TRUE,sep=',')
spplot.ab.t.z<-spplot.ab.t.z[,-1]

#categorical data (reservation and community type):
setwd("~/Documents/Research/PCAP/PCAPdata")
cats<-read.csv("2015 categorical data.csv",header=TRUE,sep=',')
categories<-matrix(as.numeric(as.factor(as.matrix(cats[,3:4]))),ncol=2)
colnames(categories)<-c("comm","res")

#set up the correlation matrix for mean abiotic factors only
data<-spplot.ab.t.z[,1:17]
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:17, 1:17, vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

#a function to run diagnostic tests on each model 
diagnostics<- function(x,print=TRUE){
  plot(x) #using Cook's distance to check for outliers
  hist(x$residuals,breaks=3) #to check for normalcy
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
#(1)
########################################################################
#make the data frame
concatenate<-data.frame(spplot.ab.t.z,categories,spplot.patterns.f)
#remove plot 10 because it's an outlier in relationship (2), and I want to be 
#able to compare models via AICc.
concatenate<-concatenate[c(1:9,11:29),]

#step (A)
#build the means-only global model
global.model<-lm(spgamma.rich ~ mOM + mP + mKppm + mMgppm + 
                   mCappm + mpH + mCEC + mK + mMg + mCa + mC + mN + mC.N +
                   mlight + ml + mo + mr, data = concatenate) 

#dredge it, subsetting the correlation matrix
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
#subset best models by showing only models within 2 points of the lowest AICc score
models<-model1[model1$delta<2,]
#plot to make sure the best model doesn't include FALSE pairs from SMAT
plot(models)

#best lm:
#check for outliers and normalcy
modA<-lm(spgamma.rich~mCa,data=concatenate)
diagnostics(modA)
#df: 1, 26
#AICc: 182.7
#multiple-r2: 54%

#step (B)
#how much of the correlation between plots is due to community type?
modB<-lme(spgamma.rich~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 24
#AICc: 194.1
#marginal-r2: 0%
#conditional-r2: 49%

#step (C1)
global.model<-lme(spgamma.rich ~ mOM + mP + mKppm + mMgppm + 
                    mCappm + mpH + mCEC + mK + mMg + mCa + mC + mN + mC.N +
                    mlight + ml + mo + mr,random=list(comm=pdDiag(~1)),method='ML', data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

#best intercepts model
modC1<-lme(spgamma.rich~mCa,random=list(comm=pdDiag(~1)),method='ML',data=concatenate)
diagnostics(modC1)
#df: 1, 23
#AICc: 185.5
#marginal-r2: 55%
#conditional-r2: 55%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[1:17]))){
  predictor<-concatenate[,i]
  model<-lme(spgamma.rich~1,random=list(comm=pdDiag(~predictor-1)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,1:17]),score)
View(results[order(results[,2]),])

#best slopes model
modC2<-lme(spgamma.rich~1,random=list(comm=pdDiag(~mCa-1)),data=concatenate,method='ML')
diagnostics(modC2)
#df: 24
#AICc: 193.3
#marginal-r2: 0%
#conditional-r2: 45%

#step (D1)
score<-NULL
for(i in 1:length(colnames(concatenate[1:17]))){
  predictor<-concatenate[,i]
  model<-lme(spgamma.rich~mCa,random=list(comm=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}#this formulation of the random parts of the model helps control for correlation between intercepts and slopes.
results<-data.frame(colnames(concatenate[,1:17]),score)
View(results[order(results[,2]),])

#best slopes constrained by intercepts model(s)
modD1.o<-lme(spgamma.rich~mCa,random=list(comm=pdDiag(~mo)),data=concatenate,method='ML')
diagnostics(modD1.o)
#df: 1, 23
#AICc: 185.6
#marginal-r2: 48%
#conditional-r2: 68%
modD1.n<-lme(spgamma.rich~mCa,random=list(comm=pdDiag(~mN)),data=concatenate,method='ML')
diagnostics(modD1.n)
#df: 1, 23
#AICc: 185.8
#marginal-r2: 36%
#conditional-r2: 75%
modD1.o.exp<-lme(spgamma.rich~mCa+mo,random=list(comm=pdDiag(~mo)),data=concatenate,method='ML')
diagnostics(modD1.o.exp)
#df: 2, 22
#AICc: 187.6
#marginal-r2: 51%
#conditional-r2: 64%
modD1.n.exp<-lme(spgamma.rich~mCa+mN,random=list(comm=pdDiag(~mN)),data=concatenate,method='ML')
diagnostics(modD1.n.exp)
#df: 2, 22
#AICc: 188.9
#marginal-r2: 34%
#conditional-r2: 73%

#step (D2)
global.model<-lme(spgamma.rich ~ mOM + mP + mKppm + mMgppm + 
                    mCappm + mpH + mCEC + mK + mMg + mCa + mC + mN + mC.N +
                    mlight + ml + mo + mr,random=list(comm=pdDiag(~mCa)),method='ML', data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

#best intercepts constrained by slopes model:
modD2<-lme(spgamma.rich~mCa,random=list(comm=pdDiag(~mCa)),method='ML',data=concatenate)
diagnostics(modD2)
#df: 1, 23
#AICc: 188.5
#marginal-r2: 55%
#conditional-r2: 55%
########################################################################
#plot best model(s):
########################################################################
#plot best AE models:
#modA has the lowest AICc score:
plot(spgamma.rich~mCa,data=concatenate,pch=16,
     main='Spring Forbs: the AE hypothesis',xlab='normalized mean calcium',
     ylab = 'plot-level richness (S)',
     cex = 1.2)
abline(coef(modA)[1],coef(modA)[2],lwd=3)
legend('topleft',legend=c('across plots'),lwd=3,col='black',pch=16,cex=1.2)

#modD1.n has the highest variance explained:
plot(spgamma.rich~mCa,data=concatenate,pch=NA,
     main='AE model (D1)',xlab='normalized mean calcium',
     ylab = 'plot-level richness (S)',cex = 1.7,xlim=c(-2,2),ylim=c(0,35))
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
abline(14.495368,5.166056,col='red',lwd=5)
#BM
points(concatenate[concatenate$comm==1,10],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.7)
abline(coef(modD1.n)[1,1],coef(modD1.n)[1,2],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,10],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.7)
abline(coef(modD1.n)[2,1],coef(modD1.n)[2,2],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,10],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.7)
abline(coef(modD1.n)[3,1],coef(modD1.n)[3,2],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,10],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modD1.n)[4,1],coef(modD1.n)[4,2],col=colors[4],lwd=5,lty=4)
legend("topleft",legend=c('across plots','Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,1,2,3,4),pch=c(NA,as.character('B'),as.character('F'),as.character('M'),
                            as.character('O')),col=c('red',colors),lwd = c(3,3,3,3,3),cex=c(1.5,1.5,1.5,1.5,1.5))

plot(spgamma.rich~mN,data=concatenate,pch=NA,
     main='AE model (D1)',xlab='normalized mean nitrogen',
     ylab = 'plot-level richness (S)',cex = 1.7,cex.axis=1.5,ylim=c(0,35))
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,12],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.7)
abline(coef(modD1.n)[1,1],coef(modD1.n)[1,3],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,12],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.7)
abline(coef(modD1.n)[2,1],coef(modD1.n)[2,3],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,12],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.7)
abline(coef(modD1.n)[3,1],coef(modD1.n)[3,3],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,12],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modD1.n)[4,1],coef(modD1.n)[4,3],col=colors[4],lwd=5,lty=4)

plot(spgamma.rich~mo,data=concatenate,pch=NA,
     main='AE model (D1)',xlab='normalized mean organic depth',
     ylab = 'plot-level richness (S)',cex = 1.2)
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,16],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1.o)[1,1],coef(modD1.o)[1,3],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,16],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1.o)[2,1],coef(modD1.o)[2,3],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,16],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1.o)[3,1],coef(modD1.o)[3,3],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,16],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1.o)[4,1],coef(modD1.o)[4,3],col=colors[4],lwd=3,lty=4)

########################################################################
#GAMMA EVENNESS
#(1)
########################################################################
#make the data frame:
concatenate<-data.frame(spplot.ab.t.z,categories,spplot.patterns.f)
#removed plot 24 due to NA values, and 10 since it was an outlier in (2)
concatenate<-concatenate[c(1:9,11:23,25:29),] 

#need richness as a variable in smat to know which abiotic factors to exclude
rich<-scale(concatenate$spgamma.rich,center=TRUE)
#data<-data.frame(rich,spplot.ab.t.z[c(1:9,11:23,25:29),1:17])
#then run smat matrix code from is.correlated

#removed mOM, mMgppm, mCappm, mpH, mMg, mCa, mC,mC.N because they are correlated with richness
data<-data.frame(spplot.ab.t.z[c(1:23,25:29),c(2:3,7:8,12,14:17)])
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat<-outer(1:dim(data)[2],1:dim(data)[2],vCorrelated,data=data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

#step (A)
global.model<-lm(log(spgamma.evens) ~ mP + mKppm + mCEC + mK + mN + 
                   mlight + ml + mo + mr, data = concatenate) 

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

modA<-lm(log(spgamma.evens)~mP+mCEC,data=concatenate)
diagnostics(modA)
#df: 2, 24
#AICc: 44.0
#adjusted-r2: 36%

#step (B)
modB<-lme(log(spgamma.evens)~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 23
#AICc: 55.8
#marginal-r2: 0%
#conditional-r2: 0%
########################################################################
#plot best model(s):
########################################################################
#plot best model modA:
vif(modA) #all <2
plot(log(spgamma.evens)~mP,data=concatenate,pch=16,col='black',cex=1.2,
     xlab='normalized mean phosphorus',ylab='species evenness (E)',
     main='Spring Forbs: evenness and phosphorus quality
     across plots',axes=FALSE)
axis(1,at=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
axis(2,at=c(-3,-2,-1.5,-1,-.5,0),labels=round(exp(c(-3,-2,-1.5,-1,-.5,0)),1))
box()
abline(coef(modA)[1],coef(modA)[2],lwd=3,col='black')
legend('bottomright',legend=c('across plots'),lwd=3,col='black',pch=16,cex=1.2)

plot(log(spgamma.evens)~mCEC,data=concatenate,pch=16,col='black',cex=1.2,
     xlab='normalized mean cation exchange capacity',ylab='species evenness (E)',
     main='Spring Forbs: evenness and cation exchange 
     capacity quality across plots',axes=FALSE)
axis(1,at=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
axis(2,at=c(-3,-2,-1.5,-1,-.5,0),labels=round(exp(c(-3,-2,-1.5,-1,-.5,0)),1))
box()
abline(coef(modA)[1],coef(modA)[3],lwd=3,col='black')
legend('topleft',legend=c('across plots'),lwd=3,col='black',pch=16,cex=1.2)


