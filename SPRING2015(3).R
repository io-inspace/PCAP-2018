#SPRING 2015
#(3) HETEROGENEITY-HETEROGENEITY RELATIONSHIP (HHR) hypothesis

#call data etc:
############################################################################
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
setwd("~/Documents/Research/PCAPdata")
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
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

#a function to run the three diagnostic tests on each model 

diagnostics<- function(x,print=TRUE){
  plot(x) #using Cook's distance to check for outliers
  hist(x$residuals,breaks=3) #to check for normalcy
  a<-coef(x) #to get model coefficientss
  b<-summary(x) #to get coefficients and stats
  c<-AICc(x)
  d<-r.squaredGLMM(x) #don't use for lms
  output<-list(a,b,c,d)
  return(output)
}

ctrl<-lmeControl(opt='optim')
options(na.action=na.fail)
###############################################################################
#BETA RICH
#(3)
###############################################################################
concatenate<-data.frame(spplot.ab.t.z,categories,spplot.patterns.f)
concatenate<-concatenate[c(1:23,25:29),] #remove NAs

#step (A)
global.model<-lm(spbeta.rich ~ vOM + vP + vKppm + vMgppm + 
                   vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN + vC.N +
                   vlight + vl + vo + vr, data = concatenate) 
model1<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model1[model1$delta<2,]
plot(models)

modA<-lm(spbeta.rich~1,data=concatenate)
diagnostics(modA)
#df: 27
#AICc: -13.9
#r2: n/a

#step (B)
modB<-lme(spbeta.rich~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 24
#AICc: -11.4
#marginal-r2: 0%
#conditional-r2: 0%
###############################################################################
#plot best model(s):
###############################################################################
plot((1-spbeta.rich)~vMg,data=concatenate,pch=16,xlab='normalized magnesium heterogeneity',ylab = 'spatial turnover of species',
     main='Spring Forbs: no support for the HHR hypothesis',cex.axis=1.5,cex=1.7,ylim=c(0,1))
legend('topleft',legend=c('plot'),pch=16,cex=1.5)
###############################################################################
#RICH HET
#(3)
###############################################################################
concatenate<-data.frame(spplot.ab.t.z,categories,spplot.patterns.f)
concatenate<-concatenate[c(1:23,25:29),] #remove NAs

#step (A)
global.model<-lm(spvar.rich ~ vOM + vP + vKppm + vMgppm + 
                   vCappm + vpH + vCEC + vK + vMg + vCa + vC + vN +
                   vC.N + vlight + vl + vo + vr, data = concatenate) 
model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)
  }))
models<-model[model$delta<2,]
plot(models)

modA<-lm(spvar.rich~1,data = concatenate)
diagnostics(modA)
#df: 27
#AICc: 244.5
#r2: n/a

#step (B)
modB<-lme(spvar.rich~1,random=list(comm=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 24
#AICc: 247.0
#marginal-r2: 0%
#conditional-r2: 4%
###############################################################################
#plot best model(s):
###############################################################################
plot(spvar.rich~vMg,data=concatenate,pch=16,xlab='normalized magnesium heterogeneity',ylab = 'richness heterogeneity',
     main='Spring Forbs: no support for the HHR hypothesis',cex=1.7,cex.axis=1.5,ylim=c(0,100))
legend('topleft',legend=c('plot'),pch=16,cex=1.5)
