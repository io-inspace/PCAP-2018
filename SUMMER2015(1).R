#SUMMER 2015
#(1) PLOT-LEVEL AVAILABLE ENERGY (AE) hypothesis

#call data etc:
########################################################################
library(MuMIn)
library(nlme)
library(lme4)
library(RColorBrewer)

setwd("~/Documents/Research/PCAP/PCAPdata/Dataframes")

#SUMMER 2015

#call raw community pattern data
summplot.patterns.a<-read.csv("summplot.patterns.allforms.csv",header=TRUE,sep=',')
summplot.patterns.a<-summplot.patterns.a[,-1]
summplot.patterns.f<-read.csv("summplot.patterns.forbs.csv",header=TRUE,sep=',')
summplot.patterns.f<-summplot.patterns.f[,-1]
summplot.patterns.o<-read.csv("summplot.patterns.otherforms.csv",header=TRUE,sep=',')
summplot.patterns.o<-summplot.patterns.o[,-1]

#call the transformed/normalized abiotic data
summplot.ab.t.z<-read.csv("summplot.ab.t.z.csv",header=TRUE,sep=',')
summplot.ab.t.z<-summplot.ab.t.z[,-1]

#call categorical data reservation and community
setwd("~/Documents/Research/PCAP/PCAPdata")
cats<-read.csv("2015 categorical data.csv",header=TRUE,sep=',')
categories<-matrix(as.numeric(as.factor(as.matrix(cats[,3:4]))),ncol=2)
colnames(categories)<-c("comm","res")

#set up the correlation matrix for mean abiotic factors only
#removed mean organic layer depth because it is extremely 0-inflated.

data<-summplot.ab.t.z[,c(1:15,17)]
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

#a function to run the three diagnostic tests on each model 
diagnostics<- function(x,print=TRUE){
  plot(x) #using Cook's distance to check for outliers
  hist(x$residuals) #to check for normalcy
  a<-coef(x) #to get model coefficients
  b<-summary(x) #to get coefficients and stats
  c<-AICc(x)
  d<-r.squaredGLMM(x)
  output<-list(a,b,c,d)
  return(output)
}

ctrl<-lmeControl(opt='optim')
options(na.action=na.fail)
########################################################################
#GAMMA RICHNESS
#(1)
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
#remove 6 so we can compare to (2) and combination models
concatenate<-concatenate[c(1:5,7:29),] 

#step (A)
global.model<-lm(summgamma.rich ~ mOM + mP + mKppm + mMgppm + 
                   mCappm + mpH + mCEC + mK + mMg + mCa + mC + mN + mC.N +
                   mlight + ml + mr, data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

#best linear model:
modA<-lm(summgamma.rich~mCa+mK+mP,data=concatenate)
diagnostics(modA)
#df: 3, 24
#AICc: 200.9
#adjusted-r2: 61%

#step (B)
modB<-lme(summgamma.rich~1,~1|as.factor(comm),data=concatenate,method='ML')
diagnostics(modB)
#df: 24
#AICc: 220.5
#marginal-r2: 0%
#conditional-r2: 30%

#step (C1)
global.model<-lme(summgamma.rich~ mOM + mP + mKppm + mMgppm + 
                    mCappm + mpH + mCEC + mK + mMg + mCa + mC + mN + mC.N +
                    mlight + ml + mr,random=list(comm=pdDiag(~1)), 
                  data = concatenate,method='ML') 

model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model1[model1$delta<2,]
plot(models)

modC1<-lme(summgamma.rich~mCa+mK+mP,random=list(comm=pdDiag(~1)), 
            data = concatenate,method='ML')
diagnostics(modC1)
#df: 3, 21
#AICc: 204.2
#marginal-r2: 66%
#conditional-r2: 66%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[1:17]))){
  predictor<-concatenate[,i]
  model<-lme(summgamma.rich~1,random=list(comm=pdDiag(~predictor-1)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}

results<-data.frame(colnames(concatenate[,1:17]),score)
View(results[order(results[,2]),])

modC2<-model<-lme(summgamma.rich~1,random=list(comm=pdDiag(~mpH-1)),data=concatenate,method='ML')
diagnostics(modC2)
#df: 24
#AICc: 215.0
#marginal-r2: 0%
#conditional-r2: 51%

#step (D1)
score<-NULL
for(i in 1:length(colnames(concatenate[1:17]))){
  predictor<-concatenate[,i]
  model<-lme(summgamma.rich~mCa+mK+mP,random=list(comm=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}

results<-data.frame(colnames(concatenate[,1:17]),score)
View(results[order(results[,2]),])

modD1<-lme(summgamma.rich~mCa+mK+mP,random=list(comm=pdDiag(~mC.N)), 
            data = concatenate,method='ML')
diagnostics(modD1)
#df: 3, 21
#AICc: 207.8
#marginal-r2: 66%
#conditional-r2: 66%
modD1.exp<-lme(summgamma.rich~mCa+mK+mP+mC.N,random=list(comm=pdDiag(~mC.N)), 
              data = concatenate,method='ML')
diagnostics(modD1.exp)
#df: 4, 20
#AICc: 211.6
#marginal-r2: 67%
#conditional-r2: 67%

#step (D2)
global.model<-lme(summgamma.rich~ mOM + mP + mKppm + mMgppm + 
                    mCappm + mpH + mCEC + mK + mMg + mCa + mC + mN + mC.N +
                    mlight + ml + mr,random=list(comm=pdDiag(~mpH)), 
                  data = concatenate,method='ML') 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model1[model1$delta<2,]
plot(models)

modD2<-lme(summgamma.rich~mCa+mK+mP,random=list(comm=pdDiag(~mpH)), 
              data = concatenate,method='ML')
diagnostics(modD2)
#df: 3, 21
#AICc: 207.8
#marginal-r2: 66%
#conditional-r2: 66%
modD2.exp<-lme(summgamma.rich~mCa+mK+mP+mpH,random=list(comm=pdDiag(~mpH)), 
              data = concatenate,method='ML')
diagnostics(modD2.exp)
#df: 4, 20
#AICc: 211.4
#marginal-r2: 67%
#conditional-r2: 67%
########################################################################
#plot best model(s)
########################################################################
plot(summgamma.rich~mCa,data=concatenate,pch=16,cex=2,col='black',
     main='AE model (A)',xlab='normalized mean calcium',
     ylab = 'plot-level richness (S)',cex.axis=1.5)
abline(coef(modA)[1],coef(modA)[2],col='black',lwd=5)
legend('topleft',legend=c('across plots'),pch=16, cex = 1.5, lwd=3)
plot(summgamma.rich~mK,data=concatenate,pch=16,col='black',cex=2,
     main='AE model (A)',xlab='normalized mean potassium',
     ylab = 'plot-level richness (S)',cex.axis=1.5)
abline(coef(modA)[1],coef(modA)[3],col='black',lwd=5)
plot(summgamma.rich~mP,data=concatenate,pch=16,col='black',cex=2,
     main='AE model (A)',xlab='normalized mean phosphorus',
     ylab = 'plot-level richness (S)',cex.axis=1.5)
abline(coef(modA)[1],coef(modA)[4],col='black',lwd=5)
########################################################################
#GAMMA EVENNESS
#(1)
########################################################################
#ok how about you remove richness and correlated terms:
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
rich<-scale(concatenate$summgamma.rich,center=TRUE)
#data<-data.frame(rich,summplot.ab.t.z[,c(1:17)])

data<-data.frame(summplot.ab.t.z[,c(2:3,7:8,14:15)])
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)
#that means richness, mOM,mMgppm, mCappm, mpH, mMg, mCa, mC, mN, mC.N, mr need to be removed.

#step (A)
global.model<-lm(summgamma.evens ~ mP + mKppm +  mCEC + mK + 
                   mlight + ml, data = concatenate) 
#+ vOM + vP + vKppm + vMgppm + 
#vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
# vlight + vl + vo + vr
model<-dredge(global.model,subset=smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

#best linear model
modA<-lm(summgamma.evens~mlight,data=concatenate)
diagnostics(modA)
#df: 1, 27
#AICc: 20.6
#multiple-r2: 13%

#step (B)
modB<-lme(summgamma.evens~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 25
#AICc: 23.5
#marginal-r2: 0%
#conditional-r2: 12%
########################################################################
#plot best model(s)
########################################################################
plot(summgamma.evens~mlight,main='Summer Forbs: mean abiotic conditions
     and species evenness',xlab='normalized mean light',
     ylab='plot-level evenness (E)',data=concatenate,pch=16,cex=1.2,ylim=c(0,1))
abline(coef(modA)[1],coef(modA)[2],lwd=3)
legend('topright',legend=c('across plots'),pch=16,lwd=1.5,cex=1.2)
########################################################################
#GAMMA EVENNESS
#uneven plots
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)

global.model<-lm(summgamma.evens[concatenate$summgamma.evens<.5] ~ mP[concatenate$summgamma.evens<.5] + mKppm[concatenate$summgamma.evens<.5] +  mCEC[concatenate$summgamma.evens<.5] + mK[concatenate$summgamma.evens<.5] + 
                   mlight[concatenate$summgamma.evens<.5] + ml[concatenate$summgamma.evens<.5], data = concatenate) 
model<-dredge(global.model,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

mod.uneven<-lm(summgamma.evens[concatenate$summgamma.evens<.5] ~ mK[concatenate$summgamma.evens<.5],data = concatenate) 
diagnostics(mod.uneven)
#df: 1, 17
#AICc: -52.2
#multiple-r2: 19%

plot(summgamma.evens[concatenate$summgamma.evens<.5]~mK[concatenate$summgamma.evens<.5],
     data=concatenate,main='Summer Forbs: mean abiotic conditions 
     and species evenness < .5',xlab='normalized mean potassium',ylab='plot-level evenness',
     pch=16,cex=1.2,ylim=c(0,1))
abline(coef(mod.uneven)[1],coef(mod.uneven)[2],lwd=3)
legend('topright',legend=c('across plots'),pch=16,cex=1.2,lwd=1.5)
########################################################################
#GAMMA EVENNESS
#even plots
########################################################################
global.model<-lm(summgamma.evens[concatenate$summgamma.evens>.5] ~ mP[concatenate$summgamma.evens>.5] + mKppm[concatenate$summgamma.evens>.5] +  mCEC[concatenate$summgamma.evens>.5] + mK[concatenate$summgamma.evens>.5] + 
                   mlight[concatenate$summgamma.evens>.5] + ml[concatenate$summgamma.evens>.5], data = concatenate) 
model<-dredge(global.model,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

mod.even<-lm(summgamma.evens[concatenate$summgamma.evens>.5] ~ mlight[concatenate$summgamma.evens>.5], 
           data = concatenate) 
diagnostics(mod.even)
#df: 1, 8
#AICc: -10.0
#multiple-r2: 58%

plot(summgamma.evens[concatenate$summgamma.evens>.5]~mlight[concatenate$summgamma.evens>.5],
     main='Summer Forbs: mean abiotic conditions and 
     species evenness > .5',xlab='normalized mean light',ylab='plot-level evenness',
     data=concatenate,pch=16,cex=1.2,ylim = c(0,1))
abline(coef(mod.even)[1],coef(mod.even)[2],lwd=3)
legend('bottomright',legend=c('across plots'),pch=16,cex=1.2,lwd=1.5)

