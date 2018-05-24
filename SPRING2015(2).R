#SPRING 2015
#(2) the HETEROGENEITY-DIVERSITY RELATIONSHIP (HDR) hypothesis

#call data etc:
###############################################################################
#libraries+
library(MuMIn)
library(nlme)
library(lme4)
library(RColorBrewer)
library(AICcmodavg)
library(car)

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

data<-spplot.ab.t.z[,18:34]
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:17, 1:17, vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

#a function to run the three diagnostic tests on each model 

diagnostics<- function(x,print=TRUE){
  plot(x) #using Cook's distance to check for outliers
  hist(x$residuals,breaks=3) #to check for normalcy
  a<-coef(x) #to get model coefficients
  b<-summary(x) #to get coefficients and stats
  c<-AICc(x)
  d<-r.squaredGLMM(x) #don't use for lm models
  output<-list(a,b,c,d)
  return(output)
}

options(na.action=na.fail)
ctrl<-lmeControl(opt='optim')
###############################################################################
#GAMMA RICH
#(2)
###############################################################################
concatenate<-data.frame(spplot.ab.t.z,categories,spplot.patterns.f)
#remove plot 10 because it's an outlier in one of the models I'd like to compare this to
concatenate<-concatenate[c(1:9,11:29),] 

#step (A)
#build the variances-only global model
global.model<-lm(spgamma.rich ~ vOM + vP + vKppm + vMgppm + 
                   vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
                   vlight + vl + vo + vr, data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

#best lm
modA<-lm(spgamma.rich~vCappm+vK,data=concatenate)
diagnostics(modA)
#df: 2, 25
#AICc: 187.9
#adjusted-r2: 45%

#step (B)
modB<-lme(spgamma.rich~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#same as in relationship 1. 

#step (C1)
global.model<-lme(spgamma.rich~ vOM + vP + vKppm + vMgppm + 
                    vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + 
                    vC.N + vlight + vl + vo + vr,random=list(comm=pdDiag(~1)), 
                  data = concatenate,method='ML') 

model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model1[model1$delta<2,]
plot(models)

#best intercepts model:
modC1<-lme(spgamma.rich~vP+vK,random=list(comm=pdDiag(~1)),data=concatenate,method='ML',control=ctrl)
diagnostics(modC1)
#df: 2, 22
#AICc: 183.4
#marginal-r2: 27% 
#conditional-r2: 69%

#step (C2)
score<-NULL
for(i in 1:length(concatenate[,1:17])){
  predictor<-concatenate[,i+17]
  model<-lme(spgamma.rich~1,random=list(comm=pdDiag(~predictor-1)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,18:34]),score)
View(results[order(results[,2]),])

#best slopes model:
modC2<-lme(spgamma.rich~1,random=list(comm=pdDiag(~vCa-1)),data=concatenate,method='ML')
diagnostics(modC2)
#df: 24
#AICc: 197.6
#marginal-r2: 0%
#conditional-r2: 57%

#step (D1)
score<-NULL
for(i in 1:length(concatenate[,1:17])){
  predictor<-concatenate[,i+17]
  model<-lme(spgamma.rich~vP+vK,random=list(comm=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,18:34]),score)
View(results[order(results[,2]),])

#best model constrained by intercepts
modD1<-lme(spgamma.rich~vP+vK,random=list(comm=pdDiag(~vCEC)),data=concatenate,method='ML')
diagnostics(modD1)
#df: 2, 22
#AICc: 183.4
#marginal-r2: 24%
#conditional-r2: 75%
modD1.exp<-lme(spgamma.rich~vP+vK+vCEC,random=list(comm=pdDiag(~vCEC)),data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 3, 21
#AICc: 186.8
#marginal-r2: 28%
#conditional-r2: 77%

#step (D2)
global.model<-lme(spgamma.rich~ vOM + vP + vKppm + vMgppm + 
                    vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + 
                    vC.N + vlight + vl + vo + vr,random=list(comm=pdDiag(~vCa)), 
                  data = concatenate,method='ML') 

model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model1[model1$delta<2,]
plot(models)

#best model constrained by slopes
modD2<-lme(spgamma.rich~vK+vP,random=list(comm=pdDiag(~vCa)),data=concatenate,method='ML',control=ctrl)
diagnostics(modD2)
#df: 2, 22
#AICc: 186.5
#marginal-r2: 26%
#conditional-r2: 73%
modD2.exp<-lme(spgamma.rich~vK+vP+vCa,random=list(comm=pdDiag(~vCa)),data=concatenate,method='ML',control=ctrl)
diagnostics(modD2.exp)
#df: 3, 21
#AICc: 190.0
#marginal-r2: 25%
#conditional-r2: 73%
###############################################################################
#plot best model(s)
###############################################################################
#plot best HDR model modD1:
vif(modD1) #all <2
plot(spgamma.rich~vP,data=concatenate,pch=NA,cex=1.2,
     main='Spring Forbs: the HDR hypothesis',
     xlab='normalized phosphorus heterogeneity',
     ylab = 'plot-level richness (S)')
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
abline(14.343047,-2.234211,col='red',lwd=3)
#BM
points(concatenate[concatenate$comm==1,19],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1)[1,1],coef(modD1)[1,2],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,19],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1)[2,1],coef(modD1)[2,2],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,19],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1)[3,1],coef(modD1)[3,2],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,19],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1)[4,1],coef(modD1)[4,2],col=colors[4],lwd=3,lty=4)
legend("topright",legend=c('across plots','Beech-Maple','Floodplain','Mixed','Oak'),
       pch=c(NA,as.character('B'),as.character('F'),as.character('M'),
             as.character('O')),col=c('red',colors),lty=c(1,1,2,3,4),lwd=c(3,1.5,1.5,1.5,1.5))

plot(spgamma.rich~vK,data=concatenate,pch=NA,cex=1.2,
     main='Spring Forbs: the HDR hypothesis',
     xlab='normalized potassium heterogeneity',
     ylab = 'plot-level richness (S)')
abline(14.343047,-3.106163,lwd=3,col='red')
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,25],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1)[1,1],coef(modD1)[1,3],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,25],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1)[2,1],coef(modD1)[2,3],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,25],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1)[3,1],coef(modD1)[3,3],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,25],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1)[4,1],coef(modD1)[4,3],col=colors[4],lwd=3,lty=4)

plot(spgamma.rich~vCEC,data=concatenate,pch=NA,cex=1.2,
     main='Spring Forbs: the HDR hypothesis',
     xlab='normalized CEC heterogeneity',
     ylab = 'plot-level richness (S)')
colors<-c('darkgreen','dodgerblue','darkorchid4','goldenrod3')
#BM
points(concatenate[concatenate$comm==1,24],
       concatenate[concatenate$comm==1,37],pch=as.character('B'),col=colors[1],cex=1.2)
abline(coef(modD1)[1,1],coef(modD1)[1,4],col=colors[1],lwd=3,lty=1)
#FP
points(concatenate[concatenate$comm==3,24],
       concatenate[concatenate$comm==3,37],pch=as.character('F'),col=colors[2],cex=1.2)
abline(coef(modD1)[2,1],coef(modD1)[2,4],col=colors[2],lwd=3,lty=2)
#M
points(concatenate[concatenate$comm==5,24],
       concatenate[concatenate$comm==5,37],pch=as.character('M'),col=colors[3],cex=1.2)
abline(coef(modD1)[3,1],coef(modD1)[3,4],col=colors[3],lwd=3,lty=3)
#O
points(concatenate[concatenate$comm==8,24],
       concatenate[concatenate$comm==8,37],pch=as.character('O'),col=colors[4],cex=1.2)
abline(coef(modD1)[4,1],coef(modD1)[4,4],col=colors[4],lwd=3,lty=4)
###############################################################################
#GAMMA EVENS
#(2)
###############################################################################
#make the data frame:
concatenate<-data.frame(spplot.ab.t.z,categories,spplot.patterns.f)
#removed plot 24 due to NA values, and 10 since it is an outlier
concatenate<-concatenate[c(1:9,11:23,25:29),] 

#need richness as a variable in smat to know which abiotic factors to exclude
rich<-scale(concatenate$spgamma.rich,center=TRUE)
#data<-data.frame(rich,spplot.ab.t.z[c(1:9,11:23,25:29),18:34])
#then run smat matrix code from is.correlated

#removed vCappm, vCEC, vK, vMg, vCa because they are correlated with richness
data<-data.frame(spplot.ab.t.z[c(1:23,25:29),c(18:21,23,28:34)])
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
global.model<-lm(log(spgamma.evens) ~ vOM + vP + vKppm + vMgppm + 
                   vpH +  vC + vN + vC.N +
                   vlight + vl + vo + vr, data = concatenate)

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

modA<-lm(log(spgamma.evens)~vl+vP,data=concatenate)
diagnostics(modA)
#df: 2, 24
#AICc: 39.4
#adjusted-r2: 46%

#step (B)
modB<-lme(log(spgamma.evens)~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 23
#AICc: 55.8
#marginal-r2: 0%
#conditional-r2: 0%
###############################################################################
#plot best model(s):
###############################################################################
#plot best model modA:
vif(modA) #all <2
plot(log(spgamma.evens)~vP,data=concatenate,pch=16,col='black',cex=1.7,
     xlab='normalized phosphorus heterogeneity',ylab='species evenness (E)',
     main='Spring Forbs: evenness and phosphorus heterogeneity 
     across plots',
     axes=FALSE)
axis(1,at=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3),cex.axis=1.5)
axis(2,at=c(-3,-2,-1.5,-1,-.5,0),labels=round(exp(c(-3,-2,-1.5,-1,-.5,0)),1),cex.axis=1.5)
box()
abline(coef(modA)[1],coef(modA)[3],lwd=5,col='black')
legend('topleft',legend=c('across plots'),lwd=3,col='black',pch=16,cex=1.5)

plot(log(spgamma.evens)~vl,data=concatenate,pch=16,col='black',cex=1.7,
     xlab='normalized litter depth heterogeneity',ylab='species evenness (E)',
     main='Spring Forbs: evenness and litter depth heterogeneity 
     across plots',axes=FALSE)
axis(1,at=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3),cex.axis=1.5)
axis(2,at=c(-3,-2,-1.5,-1,-.5,0),labels=round(exp(c(-3,-2,-1.5,-1,-.5,0)),1),cex.axis=1.5)
box()
abline(coef(modA)[1],coef(modA)[2],lwd=3,col='black')
legend('topright',legend=c('across plots'),lwd=3,col='black',pch=16,cex=1.5)
