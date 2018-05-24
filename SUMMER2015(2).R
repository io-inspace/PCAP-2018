#SUMMER 2015
#(2) the HETEROGENEITY-DIVERSITY RELATIONSHIP (HDR) hypothesis

#call data etc:
########################################################################
library(MuMIn)
library(nlme)
library(lme4)
library(RColorBrewer)

setwd("~/Documents/Research/PCAP/PCAPbckup7142017/Dataframes")
summplot.patterns.a<-read.csv("summplot.patterns.allforms.csv",
                              header=TRUE,sep=",")
summplot.patterns.a<-summplot.patterns.a[,-1]
summplot.patterns.f<-read.csv("summplot.patterns.forbs.csv",header=TRUE,
                              sep=",")
summplot.patterns.f<-summplot.patterns.f[,-1]

summplot.ab.t.z<-read.csv("summplot.ab.t.z.csv",header=TRUE,sep=',')
summplot.ab.t.z<-summplot.ab.t.z[,-1]

#call categorical data reservation and community
setwd("~/Documents/Research/PCAP")
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
#(2)
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
concatenate<-concatenate[c(1:5,7:29),]

#step (A)
global.model<-lm(summgamma.rich ~ vOM + vP + vKppm + vMgppm + 
                   vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
                   vlight + vl + vr, data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

modA<-lm(summgamma.rich~vpH + vMg,data=concatenate)
diagnostics(modA)
#df: 2, 25
#AICc: 215.6
#adjusted-r2: 31%

#step (B)
modB<-lme(summgamma.rich~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 24
#AICc: 220.5
#marginal-r2: 0%
#conditional-r2: 30%

#step (C1)
global.model<-lme(summgamma.rich~ vOM + vP + vKppm + vMgppm + 
                    vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + 
                    vC.N + vlight + vl + vr,random=list(comm=pdDiag(~1)), 
                  data = concatenate,method='ML') 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model1[model1$delta<5,]
plot(models)

modC1<-lme(summgamma.rich~vN+vpH,random=list(comm=pdDiag(~1)),data=concatenate,
            method='ML',control=ctrl)
diagnostics(modC1)
#df: 2, 22
#AICc: 208.9
#marginal-r2: 28%
#conditional-r2: 74%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[1:17]))){
  predictor<-concatenate[,i+17]
  model<-lme(summgamma.rich~1,random=list(comm=pdDiag(~predictor-1)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,18:34]),score)
View(results[order(results[,2]),])

modC2<-lme(summgamma.rich~1,random=list(comm=pdDiag(~vP-1)),data=concatenate,method='ML')
diagnostics(modC2)
#df: 24
#AICc: 224.9
#marginal-r2: 0%
#conditional-r2: 9%

#step (D1)
score<-NULL
for(i in 1:length(colnames(concatenate[1:17]))){
  predictor<-concatenate[,i+17]
  model<-lme(summgamma.rich~vpH+vN,random=list(comm=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}

results<-data.frame(colnames(concatenate[,18:34]),score)
View(results[order(results[,2]),])

modD1<-lme(summgamma.rich~vpH+vN,random=list(comm=pdDiag(~vOM)),data=concatenate,method='ML')
diagnostics(modD1)
#df: 2, 22
#AICc: 210.5
#marginal-r2: 32%
#conditional-r2: 81%
modD1.exp<-lme(summgamma.rich~vpH+vN+vOM,random=list(comm=pdDiag(~vOM)),data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 3, 21
#AICc: 213.9
#marginal-r2: 32%
#conditional-r2: 81%

#step (D2)
global.model<-lme(summgamma.rich~ vOM + vP + vKppm + vMgppm + 
                    vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + 
                    vC.N + vlight + vl + vr,random=list(comm=pdDiag(~vP)), 
                  data = concatenate,method='ML') 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model1[model1$delta<5,]
plot(models)

modD2<-lme(summgamma.rich~vN+vpH,random=list(comm=pdDiag(~vP)),data=concatenate,
              method='ML',control=ctrl)
diagnostics(modD2)
#df: 2, 22
#AICc: 211.6
#marginal-r2: 28%
#conditional-r2: 75%
modD2.exp<-lme(summgamma.rich~vN+vpH+vP,random=list(comm=pdDiag(~vP)),data=concatenate,
              method='ML',control=ctrl)
diagnostics(modD2.exp)
#df: 3, 21
#AICc: 215.2
#marginal-r2: 28%
#conditional-r2: 75%
########################################################################
#plot best model(s)
########################################################################
plot(summgamma.rich~vN,data=concatenate,pch=NA,cex=1.2,
     main='HDR model (C1)',
     xlab='normalized nitrogen heterogeneity',
     ylab = 'plot-level richness (S)',cex.axis=1.5)
abline(16.638795,4.384174,lwd=3,col='red')
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,29],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.7)
abline(coef(modC1)[1,1],coef(modC1)[1,2],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,29],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.7)
abline(coef(modC1)[2,1],coef(modC1)[2,2],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,29],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.7)
abline(coef(modC1)[3,1],coef(modC1)[3,2],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,29],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modC1)[4,1],coef(modC1)[4,2],col=colors[4],lwd=5,lty=4)

plot(summgamma.rich~vpH,data=concatenate,pch=NA,cex=1.7,
     main='HDR model (C1)',
     xlab='normalized pH heterogeneity',
     ylab = 'plot-level richness (S)',cex.axis=1.5)
abline(16.638795,6.049479,lwd=5,col='red')
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,23],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.7)
abline(coef(modC1)[1,1],coef(modC1)[1,3],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,23],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.7)
abline(coef(modC1)[2,1],coef(modC1)[2,3],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,23],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.7)
abline(coef(modC1)[3,1],coef(modC1)[3,3],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,23],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modC1)[4,1],coef(modC1)[4,3],col=colors[4],lwd=5,lty=4)
legend("topleft",legend=c('across plots','Beech-Maple','Floodplain','Mixed','Oak'),
       lty=c(1,1,2,3,4),pch=c(NA,as.character('B'),as.character('F'),as.character('M'),
                              as.character('O')),col=c('red',colors),lwd=c(3,3,3,3,3),cex=1.5)
########################################################################
#GAMMA EVENNESS
#(2)
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
rich<-scale(concatenate$summgamma.rich,center=TRUE)
#data<-data.frame(rich,summplot.ab.t.z[,c(18:32,34)])

data<-data.frame(summplot.ab.t.z[,c(18:25,28:32,34)])
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)
#vMg and vCa are correlated with richness

#step (A)
global.model<-lm(summgamma.evens~ vOM + vP + vKppm + vMgppm + 
                   vCappm + vpH + vCEC + vK + vC + vN + vC.N +
                   vlight + vl + vr, data = concatenate) 
model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

modA<-lm(summgamma.evens~vKppm+vl,data=concatenate)
diagnostics(modA)
#df: 2, 26
#AICc: 20.1
#adjusted-r2: 16%

modB<-lme(summgamma.evens~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 25
#AICc: 23.5
#marginal-r2: 0%
#conditional-r2: 11%
########################################################################
#plot best model(s)
########################################################################
plot(summgamma.evens~vKppm,data=concatenate,pch=16,cex=1.2,
     main='Summer Forbs: abiotic heterogeneity and 
     species evenness',
     xlab='normalized potassium (ppm) heterogeneity',
     ylab='plot-level evenness (E)',ylim=c(0,1))
abline(coef(modA)[1],coef(modA)[2],lwd=3)
plot(summgamma.evens~vl,data=concatenate,pch=16,cex=1.2,
     main='Summer Forbs: abiotic heterogeneity and 
     species evenness',
     xlab='normalized litter depth heterogeneity',ylab='plot-level evenness (E)',ylim=c(0,1))
abline(coef(modA)[1],coef(modA)[3],lwd=3)
legend('topright',legend=c('across plots'),pch=16,cex=1.2,lwd=1.5)
########################################################################
#GAMMA EVENNESS
#uneven plots
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
concatenate<-concatenate[concatenate$summgamma.evens<.5,]

global.model<-lm(summgamma.evens ~ vOM + vP + vKppm + vMgppm + 
                   vCappm + vpH + vCEC + vK + vC + vN + vC.N +
                   vlight + vl + vr, data = concatenate) 
model<-dredge(global.model,subset=smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

mod.uneven<-lm(summgamma.evens ~ vl+vN, data = concatenate)
diagnostics(mod.uneven)
#df: 2, 16
#AICc: -52.5
#adjusted-r2: 25%

plot(summgamma.evens~vl,data=concatenate,col=1,pch=16,cex=1.2,
     main='Summer Forbs: abiotic heterogeneity
     and species evenness < .5',xlab='normalized litter depth heterogeneity',
     ylab='plot-level evenness (E)',ylim=c(0,1))
abline(coef(mod.uneven)[1],coef(mod.uneven)[2],lwd=3)
legend('topright',legend=c('across plots'),col='black',pch=16,cex=1.2,lwd=1.5)
plot(summgamma.evens~vN,data=concatenate,col=1,pch=16,cex=1.2,
     main='Summer Forbs: abiotic heterogeneity
     and species evenness < .5',xlab='normalized nitrogen heterogeneity',
     ylab='plot-level evenness (E)',ylim=c(0,1))
abline(coef(mod.uneven)[1],coef(mod.uneven)[3],lwd=3)
########################################################################
#GAMMA EVENNESS
#even plots
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
concatenate<-concatenate[concatenate$summgamma.evens>.5,]

global.model<-lm(summgamma.evens ~ vOM + vP + vKppm + vMgppm + 
                   vCappm + vpH + vCEC + vK + vC + vN + vC.N +
                   vlight + vl + vr, data = concatenate) 
model<-dredge(global.model,subset=smat,m.lim=c(1,3),extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

mod.even<-lm(summgamma.evens ~ vl+vlight, data = concatenate)
diagnostics(mod.even)
#df: 2, 7
#AICc: -12.3
#adjusted-r2: 76%

plot(summgamma.evens~vl,data=concatenate,col=1,pch=16,cex=1.7,
     main='Summer Forbs: abiotic heterogeneity
     and species evenness > .5',xlab='normalized litter depth heterogeneity',
     ylab='plot-level evenness (E)',ylim=c(0,1),cex.axis=1.5)
abline(coef(mod.even)[1],coef(mod.even)[2],lwd=5)
legend('bottomleft',legend=c('across even plots'),pch=16,cex=1.5,lwd=3)
plot(summgamma.evens~vlight,data=concatenate,col=1,pch=16,cex=1.7,
     main='Summer Forbs: abiotic heterogeneity
     and species evenness > .5',xlab='normalized light heterogeneity',
     ylab='plot-level evenness (E)',ylim=c(0,1),cex.axis=1.5)
abline(coef(mod.even)[1],coef(mod.even)[3],lwd=5)
legend('bottomleft',legend=c('across even plots'),pch=16,cex=1.5,lwd=3)
