#SUMMER 2015
#(4) subplot-level AVAILABLE ENERGY (AE) hypothesis

#call data etc.
#######################################################################
library(nlme)
library(MuMIn)
library(lme4)
library(RColorBrewer)

setwd("~/Documents/Research/PCAP/PCAPdata")
categories<-read.csv("2015 categorical data.csv",header=TRUE,sep=',')
categories<-matrix(as.numeric(as.factor(as.matrix(categories[,3:4]))),ncol=2)
colnames(categories)<-c("comm","res")
comm<-as.matrix(rep(categories[,1],each=4))
comm<-comm[c(1:22,25:116)]
res<-as.matrix(rep(categories[,2],each=4))
res<-res[c(1:22,25:116),]

setwd("~/Documents/Research/PCAP/PCAPdata/Dataframes")
summalpha.patterns.a<-read.csv("summalpha.patterns.allforms.csv",
                               header=TRUE,sep=',')
summalpha.patterns.a<-summalpha.patterns.a[,-1]
summalpha.patterns.a<-data.frame(summalpha.patterns.a,comm,res)
summalpha.patterns.f<-read.csv("summalpha.patterns.forbs.csv",
                               header=TRUE,sep=',')
summalpha.patterns.f<-summalpha.patterns.f[,-1]
summalpha.patterns.f<-data.frame(summalpha.patterns.f,comm,res)
#NAs in 27, 91:92,94,103 & 111

summalpha.ab<-read.csv("summalpha.ab.t.z.csv",header=TRUE,sep=',')
summalpha.ab<-summalpha.ab[,c(2,4:20)]

#set up the correlation matrix for mean abiotic factors only
data<-summalpha.ab[,c(2:16,18)]
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:16, 1:16, vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

diagnostics<- function(x,print=TRUE){
  plot(x) #using Cook's distance to check for outliers
  hist(x$residuals,breaks=3) #to check for normalcy
  a<-coef(x) #to get model weights
  b<-summary(x) #to get coefficients and stats
  c<-AICc(x)
  d<-r.squaredGLMM(x)
  output<-list(a,b,c,d)
  return(output)
}

ctrl<-lmeControl(opt='optim')
options(na.action=na.fail)
#######################################################################
#ALPHA RICHNESS
#(4)
#######################################################################
concatenate<-data.frame(summalpha.ab,summalpha.patterns.f)
#######################################################################
#modeling subplots in plots
#######################################################################
#step (A)
global.model<-lm(sqrt(summalpha.rich) ~ OM + P + Kppm + Mgppm + Cappm + 
                    pH + CEC + K + Mg + Ca + C + N + C.N + 
                    summalpha.light + summalpha.l + 
                    summalpha.r,data=concatenate)
model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modA<-lm(sqrt(summalpha.rich)~Ca+K+OM+P,data=concatenate)
diagnostics(modA)
#df: 4, 109
#AICc: 310.5
#adjusted-r2: 62%

#step (B)
modB<-lme(sqrt(summalpha.rich)~1,random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 85
#AICc: 271.5
#marginal-r2: 0%
#conditional-r2: 90%

#step (C1)
global.model<-lme(sqrt(summalpha.rich) ~ OM + P + Kppm + Mgppm + Cappm + 
                    pH + CEC + K + Mg + Ca + C + N + C.N + 
                    summalpha.light + summalpha.l + 
                    summalpha.r,random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')
model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modC1<-lme(sqrt(summalpha.rich)~Ca+K+OM,random=list(plot.index=pdDiag(~1)),
        data=concatenate,method='ML')
diagnostics(modC1)
#df: 3, 82
#AICc: 249.8
#marginal-r2: 37%
#conditional-r2: 85%

#mod (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(sqrt(summalpha.rich)~1,random=list(plot.index=pdDiag(~predictor-1)),data=concatenate,control='ctrl',method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

#disregarding: organic layer depth
modC2<-lme(sqrt(summalpha.rich)~1,random=list(plot.index=pdDiag(~pH-1)),data=concatenate,control='ctrl',method='ML')
diagnostics(modC2)
#df: 85
#AICc: 321.4
#marginal-r2: 0%
#conditional-r2: 84%

#step (D1)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(sqrt(summalpha.rich)~Ca+K+OM,random=list(plot.index=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}

results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

modD1<-lme(sqrt(summalpha.rich)~Ca+K+OM,random=list(plot.index=pdDiag(~K)),data=concatenate,method='ML')
diagnostics(modD1)
#df: 3, 82
#AICc: 246.8
#marginal-r2: 41%
#conditional-r2: 88%

#step (D2)
global.model<-lme(sqrt(summalpha.rich) ~ OM + P + Kppm + Mgppm + Cappm + 
                    pH + CEC + K + Mg + Ca + C + N + C.N + 
                    summalpha.light + summalpha.l + 
                    summalpha.r,random=list(plot.index=pdDiag(~pH)),data=concatenate,method='ML')
model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modD2<-lme(sqrt(summalpha.rich)~Ca+K+OM,random=list(plot.index=pdDiag(~pH)),
              data=concatenate,control='ctrl',method='ML')
diagnostics(modD2)
#df: 3, 82
#AICc: 252.1
#marginal-r2: 37%
#conditional-r2: 85%
modD2.exp<-lme(sqrt(summalpha.rich)~Ca+K+OM+pH,random=list(plot.index=pdDiag(~CEC)),
              data=concatenate,method='ML')
diagnostics(modD2.exp)
#vif is large (>2)
#df: 4, 81
#AICc: 248.7
#marginal-r2: 42%
#conditional-r2: 88%
#######################################################################
#plot best subplot in plot model(s):
#######################################################################
#best model is modD1
plot(sqrt(summalpha.rich)~Ca,data=concatenate,pch=NA,
     main='plot model D1',
     xlab='normalized calcium',axes=FALSE)
text(concatenate$Ca,sqrt(concatenate$summalpha.rich),
     labels=concatenate$plot.index,cex=1)
axis(1,at=c(-1.5,-1,-.5,0,.5,1,1.5),labels=c(-1.5,-1,-.5,0,.5,1,1.5),cex.axis=1.5)
axis(2,at=seq(0,7,1),labels=round(seq(0,7,1)^2,1),cex.axis=1.5)
box()
abline(2.8191456,0.7300710,col='red',lwd=5)
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,2])}
legend('topleft',legend=c('across subplots','within plots'),col=c('red','black'),cex=1.5,lwd=c(3,1),pch=c(NA,NA))

plot(sqrt(summalpha.rich)~OM,data=concatenate,pch=NA,
     main='plot model D1',ylab='subplot-level richness',
     xlab='normalized % organic matter',axes=FALSE)
text(concatenate$OM,sqrt(concatenate$summalpha.rich),
     labels=concatenate$plot.index,cex=1)
axis(1,at=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3),cex.axis=1.5)
axis(2,at=seq(0,7,1),labels=round(seq(0,7,1)^2,1),cex.axis=1.5)
box()
abline(2.8191456,-0.2681773,col='red',lwd=5)
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,4])}
legend('topright',legend=c('across plots','across subplots'),col=c(1,2),lwd=c(1,3),pch=c(NA,NA),cex=1.5)

plot(sqrt(summalpha.rich)~K,data=concatenate,pch=NA,
     main='plot model D1',ylab='subplot-level richness',
     xlab='normalized potassium',axes=FALSE)
text(concatenate$K,sqrt(concatenate$summalpha.rich),labels=concatenate$plot.index,cex=1)
abline(2.8191456,-0.2568621,col='red',lwd=5)
axis(1,at=c(-2,-1,0,1,2),labels=c(-2,-1,0,1,2),cex.axis=1.5)
axis(2,at=seq(0,7,1),labels=round(seq(0,7,1)^2,1),cex.axis=1.5)
box()
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,3])}
legend('topright',legend=c('across subplots','within plots'),col=c('red','black'),lwd=c(3,1),pch=c(NA,NA),cex=1.5)
#######################################################################
#modeling subplots in plots nested in forests
#######################################################################
#step (B)
modB<-lme(sqrt(summalpha.rich)~1,~1|comm/plot.index,data=concatenate,method='ML')
diagnostics(modB)
#df: 85
#AICc: 264.6
#marginal-r2: 0%
#conditional-r2: 89%

#step (C1)
global.model<-lme(sqrt(summalpha.rich) ~ OM + P + Kppm + Mgppm + Cappm + 
                    pH + CEC + K + Mg + Ca + C + N + C.N + 
                    summalpha.light + summalpha.l + 
                    summalpha.r, ~1|comm/plot.index,data=concatenate,method='ML')
model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modC1<-lme(sqrt(summalpha.rich)~Ca+K+OM,~1|comm/plot.index,data=concatenate,method='ML')
diagnostics(modC1)
#df: 3, 82
#AICc: 249.8
#marginal-r2: 28%  
#conditional-r2: 85%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(sqrt(summalpha.rich)~1,~predictor-1|comm/plot.index,data=concatenate,control=ctrl,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

#disregarding organic layer depth:
modC2<-lme(sqrt(summalpha.rich)~1,~pH-1|comm/plot.index,data=concatenate,control=ctrl,method='ML')
diagnostics(modC2)
#df: 85
#AICc: 311.2
#marginal-r2: 0%
#conditional-r2: 83%

#step (D1)

modD1<-lme(sqrt(summalpha.rich)~Ca+K+OM,~OM|comm/plot.index,data=concatenate,control=ctrl,method='ML')
AICc(modD1)

#score<-NULL
#for(i in 1:length(colnames(concatenate[2:18]))){
#  predictor<-concatenate[,i+1]
#  model<-lme(sqrt(summalpha.rich)~Ca+K+OM,~predictor|comm/plot.index,data=concatenate,control=ctrl,method='ML')
#  score[i]<-AICc(model)
#}
#results<-data.frame(colnames(concatenate[,2:18]),score)
#View(results[order(results[,2]),])

modD1<-lme(sqrt(summalpha.rich)~Ca+K+OM,~Cappm|comm/plot.index,data=concatenate,control=ctrl,method='ML')
diagnostics(modD1)
#df: 3, 82
#AICc: 247.3
#marginal-r2: 37%
#conditional-r2: 88%
modD1.exp<-lme(sqrt(summalpha.rich)~Ca+K+OM+Cappm,~Cappm|comm/plot.index,data=concatenate,control=ctrl,method='ML')
diagnostics(modD1.exp)
#vif is large
#df: 4, 81
#AICc: 249.8
#marginal-r2: 37%
#conditional-r2: 87%

#step (D2)
#do this iteratively
model<-lme(sqrt(summalpha.rich) ~ OM, ~pH|comm/plot.index,control=ctrl,data=concatenate,method='ML')


#global.model<-lme(sqrt(summalpha.rich) ~ OM + P + Kppm + Mgppm + Cappm + 
#                    pH + CEC + K + Mg + Ca + C + N + C.N + 
#                    summalpha.light + summalpha.l + 
#                    summalpha.r, ~pH|comm/plot.index,data=concatenate,method='ML')
#model<-dredge(global.model,subset = smat,extra = list(
#  "R^2", "*" = function(x) {
#    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
#  }))
#models<-model[model$delta<2,]
#plot(models)

modD2<-lme(sqrt(summalpha.rich) ~ summalpha.l, ~pH|comm/plot.index,data=concatenate,method='ML',control=ctrl)
diagnostics(modD2)
#df: 1, 84
#AICc: 259.8
#marginal-r2: 1%
#conditional-r2: 90%
modD2.exp<-lme(sqrt(summalpha.rich) ~ summalpha.l+pH, ~pH|comm/plot.index,data=concatenate,method='ML',control=ctrl)
diagnostics(modD2.exp)
#df: 2, 83
#AICc: 259.2
#marginal-r2: 14%
#conditional-r2: 87%

#######################################################################
#plot best subplot in plot nested in forests model(s):
#######################################################################
plot(sqrt(summalpha.rich)~Ca,data=concatenate,pch=NA,
     main='Summer Forbs: the subplot-level AE hypothesis
     across subplots and forests',ylab='subplot-level richness',
     xlab='normalized calcium',axes=FALSE)
text(concatenate$Ca,sqrt(concatenate$summalpha.rich),
     labels=concatenate$plot.index,cex=.7)
axis(1,at=c(-1.5,-1,-.5,0,.5,1,1.5),labels=c(-1.5,-1,-.5,0,.5,1,1.5))
axis(2,at=seq(0,7,1),labels=round(seq(0,7,1)^2,1))
box()
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,2])}
abline(2.5888022 ,0.6681532 ,col='red',lwd=3)
legend("topleft",legend=c('across subplots','within forests'),lwd=c(3,1.5),col=c('red','black'))

plot(sqrt(summalpha.rich)~K,data=concatenate,pch=NA,
     main='Summer Forbs: the subplot-level AE hypothesis
     across subplots and forests',ylab='subplot-level richness',
     xlab='normalized potassium',axes=FALSE)
text(concatenate$K,sqrt(concatenate$summalpha.rich),
     labels=concatenate$plot.index,cex=.7)
axis(1,at=c(-2,-1,0,1,2),labels=c(-2,-1,0,1,2))
axis(2,at=seq(0,7,1),labels=round(seq(0,7,1)^2,1))
box()
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,3])}
abline(2.5888022 ,-0.2711517,col='red',lwd=3)
legend("topright",legend=c('across subplots','within forests'),lwd=c(3,1.5),col=c('red','black'))

plot(sqrt(summalpha.rich)~OM,data=concatenate,pch=NA,
     main='Summer Forbs: the subplot-level AE hypothesis
     across subplots and forests',ylab='subplot-level richness',
     xlab='normalized % organic matter',axes=FALSE)
text(concatenate$OM,sqrt(concatenate$summalpha.rich),
     labels=concatenate$plot.index,cex=.7)
axis(1,at=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3))
axis(2,at=seq(0,7,1),labels=round(seq(0,7,1)^2,1))
box()
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,4])}
abline(2.5888022 ,-0.2648193,col='red',lwd=3)
legend("topright",legend=c('across subplots','within forests'),lwd=c(3,1.5),col=c('red','black'))

plot(sqrt(summalpha.rich)~Cappm,data=concatenate,pch=NA,
     main='Summer Forbs: the subplot-level AE hypothesis
     across subplots and forests',ylab='subplot-level richness',
     xlab='normalized calcium (ppm)',axes=FALSE)
text(concatenate$Cappm,sqrt(concatenate$summalpha.rich),
     labels=concatenate$plot.index,cex=.7)
axis(1,at=c(-3,-2,-1,0,1,2),labels=c(-3,-2,-1,0,1,2))
axis(2,at=seq(0,7,1),labels=round(seq(0,7,1)^2,1))
box()
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,5])}
legend("topleft",legend=c('across subplots','within forests'),lwd=c(3,1.5),col=c('red','black'))
#######################################################################
#ALPHA EVENNESS
#(4)
#######################################################################
concatenate<-data.frame(summalpha.ab,summalpha.patterns.f)
#NAs in 27, 91:92,94,103 & 111
concatenate<-concatenate[c(1:26,28:90,93,95:102,104:110,112:114),]
#rich<-scale(concatenate$summalpha.rich,center=TRUE)
#data<-data.frame(rich,summalpha.ab[c(1:26,28:90,93,95:102,104:110,112:114),2:18])

data<-summalpha.ab[c(1:26,28:90,93,95:102,104:110,112:114),c(2:4,8:9,12:16,18)]
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)
#don't care about Mgppm, Cappm, pH, Mg, Ca, because they're correlated with richness
#######################################################################
#modeling subplots in plots
#######################################################################
#step (A)
global.model<-lm(summalpha.evens ~ OM + P + Kppm + K +CEC + C + N + 
                    C.N + 
                    summalpha.light + summalpha.l + 
                    summalpha.r,data=concatenate)

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modA<-lm(summalpha.evens~C.N+CEC+summalpha.light+summalpha.l,data=concatenate)
diagnostics(modA)
#df: 4, 103
#AICc: 86.0
#adjusted-r2: 22%

#step (B)
modB<-lme(summalpha.evens~1,random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 79
#AICc: 64.3
#marginal-r2: 0%
#conditional-r2: 61%

#step (C1)
global.model<-lme(summalpha.evens ~ OM + P + Kppm + K +CEC + C + N + 
                    C.N + 
                    summalpha.light + summalpha.l + 
                    summalpha.r,random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modC1<-lme(summalpha.evens~C.N+summalpha.l+summalpha.light,
        random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modC1)
#df: 3, 76
#AICc: 61.1
#marginal-r2: 12%
#conditional-r2: 56%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(summalpha.evens~1,random=list(plot.index=pdDiag(~predictor-1)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

modC2<-lme(summalpha.evens~1,random=list(plot.index=pdDiag(~N-1)),data=concatenate,method='ML')
diagnostics(modC2)
#df: 79
#AICc: 95.2
#marginal-r2: 0%
#conditional-r2: 40%

#step (D1)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(summalpha.evens~C.N+summalpha.l+summalpha.light,random=list(plot.index=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

modD1<-lme(summalpha.evens~C.N+summalpha.l+summalpha.light,
           random=list(plot.index=pdDiag(~C.N)),data=concatenate,method='ML')
diagnostics(modD1)
#df: 3, 76
#AICc: 60.4
#marginal-r2: 16%
#conditional-r2: 65%

#step (D2)
global.model<-lme(summalpha.evens ~ OM + P + Kppm + K +CEC + C + N + 
                    C.N + 
                    summalpha.light + summalpha.l + 
                    summalpha.r,random=list(plot.index=pdDiag(~N)),data=concatenate,method='ML')

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modD2<-lme(summalpha.evens~summalpha.light+C.N,
              random=list(plot.index=pdDiag(~N)),data=concatenate,method='ML')
diagnostics(modD2)
#df: 2, 77
#AICc: 61.3
#marginal-r2: 10%
#conditional-r2: 61%
modD2.exp<-lme(summalpha.evens~summalpha.light+N+C.N,random=list(plot.index=pdDiag(~N)),data=concatenate,method='ML')
diagnostics(modD2.exp)
#df: 3, 76
#AICc: 63.2
#marginal-r2: 11%
#conditional-r2: 62%
#######################################################################
#plot best subplot in plot models:
#######################################################################
plot(summalpha.evens~C.N,data=concatenate,pch=NA,
     main='Summer Forbs: the subplot-level AE hypothesis 
     across subplots and plots',ylab='subplot-level evenness',
     xlab='normalized carbon:nitrogen')
text(concatenate$C.N,concatenate$summalpha.evens,
     labels=concatenate$plot.index,cex=.7)
abline(0.6147737,0.1094161,col='red',lwd=3)
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,2])}
legend('bottomright',legend=c('across subplots','within plots'),col=c('red','black'),lwd=c(3,1),pch=c(NA,NA))

plot(summalpha.evens~summalpha.light,data=concatenate,pch=NA,
     main='Summer Forbs: the subplot-level AE hypothesis 
     across subplots and plots',ylab='subplot-level evenness',
     xlab='normalized light')
text(concatenate$summalpha.light,concatenate$summalpha.evens,
     labels=concatenate$plot.index,cex=.7)
abline(0.6147737,-0.0755398,col='red',lwd=3)
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,4])}
legend('bottomright',legend=c('across subplots','within plots'),col=c('red','black'),lwd=c(3,1),pch=c(NA,NA))

plot(summalpha.evens~summalpha.l,data=concatenate,pch=NA,
     main='Summer Forbs: the subplot-level AE hypothesis 
     across subplots and plots',ylab='subplot-level evenness',
     xlab='normalized litter depth')
text(concatenate$summalpha.l,concatenate$summalpha.evens,
     labels=concatenate$plot.index,cex=.7)
abline(0.6147737,0.0625866 ,col='red',lwd=3)
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],coef(modD1)[i,3])}
#######################################################################
#modeling subplots in plots nested in forests
#######################################################################
#step (B)
modB<-lme(summalpha.evens~1,~1|comm/plot.index,data=concatenate,method='ML')
diagnostics(modB)
#df: 79
#AICc: 55.8
#marginal-r2: 0%
#conditional-r2: 60%

#step (C1)
global.model<-lme(summalpha.evens ~ OM + P + Kppm + CEC + K + C + N + 
                    C.N + 
                    summalpha.light + summalpha.l + 
                    summalpha.r,~1|comm/plot.index,data=concatenate,method='ML')
model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modC1<-lme(summalpha.evens~summalpha.light,~1|comm/plot.index,data=concatenate,method='ML')
diagnostics(modC1)
#df: 1, 78
#AICc: 52.2
#marginal-r2: 5%
#conditional-r2: 59%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(summalpha.evens~1,~predictor-1|comm/plot.index,control=ctrl,data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

modC2<-lme(summalpha.evens~1,~C.N-1|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modC2)
#df: 79
#AICc: 80.2
#marginal-r2: 0%
#conditional-r2: 28%

#step (D1)
#do these iteratively
model<-lme(summalpha.evens~summalpha.light,~OM|comm/plot.index,control=ctrl,data=concatenate,method='ML')
AICc(model)

#score<-NULL
#for(i in 1:length(colnames(concatenate[2:18]))){
#  predictor<-concatenate[,i+1]
#  model<-lme(summalpha.evens~summalpha.light,~predictor|comm/plot.index,control=ctrl,data=concatenate,method='ML')
#  score[i]<-AICc(model)
#}
#results<-data.frame(colnames(concatenate[,2:18]),score)
#View(results[order(results[,2]),])

modD1<-lme(summalpha.evens~summalpha.light,~C.N|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modD1)
#df: 1, 78
#AICc: 57.5
#marginal-r2: 6%
#conditional-r2: 63%
modD1.exp<-lme(summalpha.evens~summalpha.light+C.N,~C.N|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 2, 77
#AICc: 59.8
#marginal-r2: 6%
#conditional-r2: 63%

#modD2
#do this iteratively
modD2<-lme(summalpha.evens ~ OM, ~C.N|comm/plot.index,control=ctrl,data=concatenate,method='ML')
AICc(modD2)
global.model<-lme(summalpha.evens ~ OM + P + Kppm + CEC + K + C + N + 
                    C.N + 
                    summalpha.light + summalpha.l + 
                    summalpha.r,~C.N|comm/plot.index,control=ctrl,data=concatenate,method='ML')
model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modD2<-lme(summalpha.evens ~ summalpha.light, ~C.N|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modD2)
#df: 1, 78
#AICc: 57.7
#marginal-r2: 6%
#conditional-r2: 63%
modD2.exp<-lme(summalpha.evens ~ summalpha.light+C.N, ~C.N|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modD2.exp)
#df: 2, 77
#AICc: 59.8
#marginal-r2: 6%
#conditional-r2: 63%

#######################################################################
#plot best subplots in plots nested in forests model:
#######################################################################
#best model is modC1
plot(summalpha.evens~summalpha.light,data=concatenate,pch=NA,
     main='Summer Forbs: the subplot-level AE hypothesis
     and species evenness across subplots and forests',ylab='subplot-level evenness',
     xlab='normalized light',cex.axis=1.5)
abline(0.6059508 ,-0.0917038  ,col='red',lwd=5)
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
text(concatenate$summalpha.light[concatenate$comm==1],concatenate$summalpha.evens[concatenate$comm==1],
     labels=concatenate$plot.index[concatenate$comm==1],cex=1,col=colors[1])
for(i in 1:6){abline(coef(modC1)[i,1],coef(modC1)[i,2],lwd=1.5,col=colors[1])}
#FP
text(concatenate$summalpha.light[concatenate$comm==3],concatenate$summalpha.evens[concatenate$comm==3],
     labels=concatenate$plot.index[concatenate$comm==3],cex=1,col=colors[2])
for(i in 7:15){abline(coef(modC1)[i,1],coef(modC1)[i,2],lwd=1.5,lty=2,col=colors[2])}
#M
text(concatenate$summalpha.light[concatenate$comm==5],concatenate$summalpha.evens[concatenate$comm==5],
     labels=concatenate$plot.index[concatenate$comm==5],cex=1,col=colors[3])
for(i in 16:20){abline(coef(modC1)[i,1],coef(modC1)[i,2],lwd=1.5,lty=3,col=colors[3])}
#OAK
text(concatenate$summalpha.light[concatenate$comm==8],concatenate$summalpha.evens[concatenate$comm==8],
     labels=concatenate$plot.index[concatenate$comm==8],cex=1,col=colors[4])
for(i in 21:29){abline(coef(modC1)[i,1],coef(modC1)[i,2],lwd=1.5,lty=4,col=colors[4])}
legend('topright',cex=1.5,legend=c('across subplots','Beech-Maple','Floodplain','Mixed','Oak'),lwd=c(3,3,3,3,3),col=c('red',colors),lty=c(1,1,2,3,4))
