#SUMMER 2015 
#combination models (AE & HDR hypotheses)

#call data etc:
########################################################################
library(MuMIn)
library(nlme)
library(lme4)
library(RColorBrewer)

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

data<-summplot.ab.t.z[,c(2,8,10,23,29)]#includes mCa mK and mP, vpH, vN
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)
#none of these are correlated with each other.

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
########################################################################
#GAMMA RICHNESS
#combination model
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
concatenate<-concatenate[c(1:5,7:29),]#because 6 is removed for the variance model.

#step (A)
global.model<-lm(summgamma.rich ~ mCa+mK+mP+vpH+vN, data = concatenate)
model1<-dredge(global.model,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<3,]
plot(models)

#best COMBINED ONE (isn't the best one, but that will already be reflected in the AICc)
modA<-lm(summgamma.rich~mCa+mK+mP+vpH,data=concatenate)
diagnostics(modA)
#df: 4, 23
#AICc: 202.7
#adjusted-r2: 62%

#step (C1)
global.model<-lme(summgamma.rich~ mK+mCa+mP+vpH+vN,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
model1<-dredge(global.model,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<4,]
plot(models)

modC1<-lme(summgamma.rich~mK+mCa+mP+vpH,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modC1)
#df: 4, 20
#AICc: 206.3
#marginal-r2: 68%
#conditional-r2: 68%

#step (C2)
score<-NULL
for(i in 1:length(concatenate[,c(1:34)])){
  predictor<-concatenate[,i]
  model<-lme(summgamma.rich~1,random=list(comm=pdDiag(~predictor-1)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,c(1:34)]),score)
View(results[order(results[,2]),])

modC2<-lme(summgamma.rich~1,random=list(comm=pdDiag(~mpH-1)),data=concatenate,method='ML')
diagnostics(modC2)
#same as in summ2015relationship1

#step (D1)
#USING THE ACTUAL BEST MODEL, which was the AE model.
score<-NULL
for(i in 1:length(concatenate[,c(18:34)])){
  predictor<-concatenate[,i+17]
  model<-lme(summgamma.rich~mK+mCa+mP,random=list(comm=pdDiag(~predictor)),data=concatenate,control='ctrl',method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,c(18:34)]),score)
View(results[order(results[,2]),])

modD1<-lme(summgamma.rich~mK+mCa+mP,random=list(comm=pdDiag(~vP)),data=concatenate,control='ctrl',method='ML')
diagnostics(modD1)
#df: 3, 21
#AICc: 206.7
#marginal-r2: 65%
#conditional-r2: 69%
modD1.exp<-lme(summgamma.rich~mK+mCa+mP+vP,random=list(comm=pdDiag(~vP)),data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 4, 20
#AICc: 210.4
#marginal-r2: 66%
#conditional-r2: 69%

#for completeness
score<-NULL
for(i in 1:length(concatenate[,c(1:34)])){
  predictor<-concatenate[,i]
  model<-lme(summgamma.rich~mK+mCa+mP+vpH,random=list(comm=pdDiag(~predictor)),data=concatenate,control='ctrl',method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,c(18:34)]),score)
View(results[order(results[,2]),])

modD1<-lme(summgamma.rich~mK+mCa+mP+vpH,random=list(comm=pdDiag(~vCEC)),data=concatenate,method='ML')
diagnostics(modD1)
#df: 4, 20
#AICc: 208.2
#marginal-r2: 59%
#conditional-r2: 81%
modD1.exp<-lme(summgamma.rich~mK+mCa+mP+vpH+vCEC,random=list(comm=pdDiag(~vCEC)),data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 5, 19
#AICc: 212.4
#marginal-r2: 59%
#conditional-r2: 81%
#since the intercepts are the same don't put different ones in the model structure:
modD1<-lme(summgamma.rich~mK+mCa+mP+vpH,random=list(comm=pdDiag(~vCEC-1)),data=concatenate,method='ML')
diagnostics(modD1)

#step (D2)
global.model<-lme(summgamma.rich~ mK+mCa+mP+vpH+vN,random=list(comm=pdDiag(~mpH)),data=concatenate,method='ML')
model1<-dredge(global.model,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<4,]
plot(models)

#not the best model, but the best combination model, since we've already got the best model listed in relationship 1.
modD2<-lme(summgamma.rich~ mK+mCa+vpH,random=list(comm=pdDiag(~mpH)),data=concatenate,method='ML')
diagnostics(modD2)
#df: 3, 21
#AICc: 210.2
#marginal-r2: 63%
#conditional-r2: 63%
modD2.exp<-lme(summgamma.rich~ mK+mCa+vpH+mpH,random=list(comm=pdDiag(~mpH)),data=concatenate,method='ML')
diagnostics(modD2.exp)
#df: 4, 20
#AICc: 212.1
#marginal-r2: 66%
#conditional-r2: 66%
########################################################################
#plot best model(s)
########################################################################
#modA
plot(summgamma.rich~mCa,data=concatenate,cex=1.2,pch=16,col=1,
     main='Summer Forbs: the AE & HDR hypothesis
     across plots',xlab='normalized mean calcium',ylab = 'plot-level richness (S)')
abline(coef(modA)[1],coef(modA)[2],lwd=3,col=1)
legend('topleft',legend=c('across plots'),cex=1.2,pch=16,lwd=1.5)
plot(summgamma.rich~mK,data=concatenate,cex=1.2,pch=16,col=1,
     main='Summer Forbs: the AE & HDR hypothesis
     across plots',xlab='normalized mean potassium',ylab = 'plot-level richness (S)')
abline(coef(modA)[1],coef(modA)[3],lwd=3,col=1)
plot(summgamma.rich~mP,data=concatenate,cex=1.2,pch=16,col=1,
     main='Summer Forbs: the AE & HDR hypothesis
     across plots',xlab='normalized mean phosphorus',ylab = 'plot-level richness (S)')
abline(coef(modA)[1],coef(modA)[4],lwd=3,col=1)
plot(summgamma.rich~vpH,data=concatenate,cex=1.2,pch=16,col=1,
     main='Summer Forbs: the AE & HDR hypothesis
     across plots',xlab='normalized pH heterogeneity',ylab = 'plot-level richness (S)')
abline(coef(modA)[1],coef(modA)[5],lwd=3,col=1)

#modD1
plot(summgamma.rich~mCa,data=concatenate,pch=NA,cex=1.2,
     main='combination model (D1)',
     xlab='normalized mean calciium',
     ylab = 'plot-level richness (S)')
abline(17.628632 ,9.567155,lwd=3,col='red')
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,10],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1)[1,1],coef(modD1)[1,3],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,10],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1)[2,1],coef(modD1)[2,3],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,10],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1)[3,1],coef(modD1)[3,3],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,10],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1)[4,1],coef(modD1)[4,3],col=colors[4],lwd=3,lty=4)

plot(summgamma.rich~mK,data=concatenate,pch=NA,cex=1.2,
     main='combination model (D1)',
     xlab='normalized mean potassium',
     ylab = 'plot-level richness (S)')
abline(17.628632 ,-4.909602,lwd=3,col='red')
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,8],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1)[1,1],coef(modD1)[1,2],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,8],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1)[2,1],coef(modD1)[2,2],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,8],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1)[3,1],coef(modD1)[3,2],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,8],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1)[4,1],coef(modD1)[4,2],col=colors[4],lwd=3,lty=4)

plot(summgamma.rich~mP,data=concatenate,pch=NA,cex=1.2,
     main='combination model (D1)',
     xlab='normalized mean phosphorus',
     ylab = 'plot-level richness (S)')
abline(17.628632 ,-2.861561,lwd=3,col='red')
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,2],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1)[1,1],coef(modD1)[1,4],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,2],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1)[2,1],coef(modD1)[2,4],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,2],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1)[3,1],coef(modD1)[3,4],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,2],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1)[4,1],coef(modD1)[4,4],col=colors[4],lwd=3,lty=4)

plot(summgamma.rich~vpH,data=concatenate,pch=NA,cex=1.2,
     main='combination model (D1)',
     xlab='normalized pH heterogeneity',
     ylab = 'plot-level richness (S)')
abline(17.628632 , 2.478387 ,lwd=3,col='red')
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,23],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1)[1,1],coef(modD1)[1,5],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,23],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1)[2,1],coef(modD1)[2,5],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,23],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1)[3,1],coef(modD1)[3,5],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,23],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1)[4,1],coef(modD1)[4,5],col=colors[4],lwd=3,lty=4)

plot(summgamma.rich~vCEC,data=concatenate,pch=NA,cex=1.7,
     main='combination model (D1)',
     xlab='normalized cation exchange capacity heterogeneity',
     ylab = 'plot-level richness (S)',cex.axis=1.5)
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,24],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.7)
abline(coef(modD1)[1,1],coef(modD1)[1,6],col=colors[1],lwd=5,lty=1)
#FP
points(concatenate[concatenate$comm==3,24],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.7)
abline(coef(modD1)[2,1],coef(modD1)[2,6],col=colors[2],lwd=5,lty=2)
#M
points(concatenate[concatenate$comm==5,24],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.7)
abline(coef(modD1)[3,1],coef(modD1)[3,6],col=colors[3],lwd=5,lty=3)
#O
points(concatenate[concatenate$comm==8,24],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.7)
abline(coef(modD1)[4,1],coef(modD1)[4,6],col=colors[4],lwd=5,lty=4)

########################################################################
#GAMMA EVENNESS
#combination models
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)

data<-data.frame(summplot.ab.t.z[,c(14,20,32)])
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

global.model<-lm(summgamma.evens~mlight+vKppm+vl,data=concatenate)

model1<-dredge(global.model,subset=smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

#step (A)
#the best one is vKppm vl, but that's all HDR.
modA<-lm(summgamma.evens~vKppm+mlight,data=concatenate)
diagnostics(modA)
#df: 2, 26
#AICc: 20.4
#adjusted-r2: 15%
########################################################################
#uneven plots
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
concatenate<-concatenate[concatenate$summgamma.evens<.5,]

global.model<-lm(summgamma.evens~mK+mlight+vKppm+vl+vN, data = concatenate) 
model<-dredge(global.model,subset=smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<5,]
plot(models)

mod.uneven<-lm(summgamma.evens~mK+vN+vl, data = concatenate) 
diagnostics(mod.uneven)
#df: 3, 15
#AICc: -56.2
#adjusted-r2: 45%
########################################################################
#even plots
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
concatenate<-concatenate[concatenate$summgamma.evens>.5,]

global.model<-lm(summgamma.evens~mK+mlight+vKppm+vl+vN, data = concatenate) 
model<-dredge(global.model,subset=smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

mod.even<-lm(summgamma.evens~mlight+vl, data = concatenate) 
diagnostics(mod.even)
#df: 2, 7
#AICc: -11.2
#adjusted-r2: 73%
########################################################################
#plot best model(s)
########################################################################
concatenate<-data.frame(summplot.ab.t.z,categories,summplot.patterns.f)
concatenate<-concatenate[concatenate$summgamma.evens<.5,]

plot(summgamma.evens~mK,data=concatenate,col=1,pch=16,cex=1.7,main=
       'Summer Forbs: species evenness across uneven plots',xlab='normalized mean potassium',
     ylab='plot-level evenness (E)',ylim=c(0,1),cex.axis=1.5)
abline(coef(mod.uneven)[1],coef(mod.uneven)[2],col='black',lwd=5)
legend('topright',legend=c('across uneven plots'),col='black',cex=1.5,pch=16,lwd=3)

plot(summgamma.evens~vN,data=concatenate,col=1,pch=16,cex=1.7,main=
       'Summer Forbs: species evenness across uneven plots',xlab='normalized nitrogen heterogeneity',
     ylab='plot-level evenness (E)',ylim=c(0,1),cex.axis=1.5)
abline(coef(mod.uneven)[1],coef(mod.uneven)[3],col=1,lwd=5)
legend('topright',legend=c('across uneven plots'),col='black',cex=1.5,pch=16,lwd=3)

plot(summgamma.evens~vl,data=concatenate,col=1,pch=16,cex=1.7,main=
       'Summer Forbs: species evenness across uneven plots',xlab='normalized litter depth heterogeneity',
     ylab='plot-level evenness (E)',ylim=c(0,1),cex.axis=1.5)
abline(coef(mod.uneven)[1],coef(mod.uneven)[4],col=1,lwd=5)
legend('topright',legend=c('across uneven plots'),col='black',cex=1.5,pch=16,lwd=3)
