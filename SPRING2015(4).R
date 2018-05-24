#SPRING 2015
#(4) the subplot-level AVAILABLE ENERGY (AE) hypothesis

#call data etc:
########################################################################
setwd("~/Documents/Research/PCAP/PCAPdata")
categories<-read.csv("2015 categorical data.csv",header=TRUE,sep=',')
categories<-matrix(as.numeric(as.factor(as.matrix(categories[,3:4]))),ncol=2)
colnames(categories)<-c("comm","res")
comm<-as.matrix(rep(categories[,1],each=4))
comm<-comm[c(1:22,25:116)]
res<-as.matrix(rep(categories[,2],each=4))
res<-res[c(1:22,25:116),]

setwd("~/Documents/Research/PCAP/PCAPdata/Dataframes")
spalpha.patterns.a<-read.csv("spalpha.patterns.allforms.csv",
                             header=TRUE,sep=',')
spalpha.patterns.a<-spalpha.patterns.a[,-1]
spalpha.patterns.a<-data.frame(spalpha.patterns.a,comm,res)
spalpha.patterns.f<-read.csv("spalpha.patterns.forbs.csv",
                             header=TRUE,sep=',')
spalpha.patterns.f<-spalpha.patterns.f[,-1]
spalpha.patterns.f<-data.frame(spalpha.patterns.f,comm,res)
#NAs in 91:94 & 113

spalpha.ab<-read.csv("spalpha.ab.t.z.csv",header=TRUE,sep=',')
spalpha.ab<-spalpha.ab[,c(2,4:20)]

#set up the correlation matrix for mean abiotic factors only
data<-spalpha.ab[,2:18]
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:17, 1:17, vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)

diagnostics<- function(x,print=TRUE){
  plot(x) #using Cook's distance to check for outliers
  hist(x$residuals) #to check for normalcy
  a<-coef(x) #to get model coefficients
  b<-summary(x) #to get coefficients and stats
  c<-AICc(x)
  d<-r.squaredGLMM(x) #do not use for lms
  output<-list(a,b,c,d)
  return(output)
}

ctrl<-lmeControl(opt='optim')
options(na.action=na.fail)
########################################################################
#ALPHA RICHNESS
#(4)
########################################################################
#modeling subplots in plots
########################################################################
concatenate<-data.frame(spalpha.ab,spalpha.patterns.f)
#step (B)
modB.plots<-lme(sqrt(spalpha.rich)~1,random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB.plots)
#df: 85
#AICc: 214.4
#marginal-r2: 0%
#conditional-r2: 90%

#step (C1)
global.model<-lme(sqrt(spalpha.rich) ~ OM + P + Kppm + Mgppm + Cappm + 
                    pH + CEC + K + Mg + Ca + C + N + C.N + 
                    spalpha.light + spalpha.l + spalpha.o + 
                    spalpha.r, random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

#best intercepts model:
modC1<-lme(sqrt(spalpha.rich)~Ca,random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modC1)
#df: 1, 84
#AICc: 202.2
#marginal-r2: 11%
#conditional-r2: 88%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(sqrt(spalpha.rich)~1,random=list(plot.index=pdDiag(~predictor-1)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

#best slopes model:
modC2<-lme(sqrt(spalpha.rich)~1,random=list(plot.index=pdDiag(~Cappm-1)),data=concatenate,method='ML')
diagnostics(modC2)
#df: 85
#AICc: 273.4
#marginal-r2: 0%
#conditional-r2: 87%

#step (D1)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(sqrt(spalpha.rich)~Ca,random=list(plot.index=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

#best model constrained by intercepts:
modD1<-lme(sqrt(spalpha.rich)~Ca,random=list(plot.index=pdDiag(~spalpha.light)),data=concatenate,method='ML')
diagnostics(modD1)
#df: 1, 84
#AICc: 204.0
#marginal-r2: 12%
#conditional-r2: 88%
modD1.exp<-lme(sqrt(spalpha.rich)~Ca+spalpha.light,random=list(plot.index=pdDiag(~spalpha.light)),data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 2, 83
#AICc: 206.0
#marginal-r2: 12%
#conditional-r2: 88%

#step (D2)
global.model<-lme(sqrt(spalpha.rich) ~ OM + P + Kppm + Mgppm + Cappm + 
                    pH + CEC + K + Mg + Ca + C + N + C.N + 
                    spalpha.light + spalpha.l + spalpha.o + 
                    spalpha.r, random=list(plot.index=pdDiag(~Cappm)),data=concatenate,method='ML')

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

#best model constrained by slopes
modD2<-lme(sqrt(spalpha.rich)~Ca,random=list(plot.index=pdDiag(~Cappm)),data=concatenate,method='ML')
diagnostics(modD2)
#df: 1, 84
#AICc: 204.3
#marginal-r2: 11%
#conditional-r2: 88%
modD2.exp<-lme(sqrt(spalpha.rich)~Ca+Cappm,random=list(plot.index=pdDiag(~Cappm)),data=concatenate,method='ML')
diagnostics(modD2.exp)
#df: 1, 84
#AICc: 205.3
#marginal-r2: 11%
#conditional-r2: 88%
########################################################################
#modeling subplots in plots nested in forests:
########################################################################
concatenate<-data.frame(spalpha.ab,spalpha.patterns.f)
#step (B)
modB.plotsinforests<-lme(sqrt(spalpha.rich)~1,~1|comm/plot.index,data=concatenate,method='ML')
diagnostics(modB.forest.plot)
#df: 85
#AICc: 211.3
#marginal-r2: 0%
#conditional-r2: 89%

#step (C1)
global.model<-lme(sqrt(spalpha.rich) ~ OM + P + Kppm + Mgppm + Cappm + 
                    pH + CEC + K + Mg + Ca + C + N + C.N + 
                    spalpha.light + spalpha.l + spalpha.o + 
                    spalpha.r, ~1|comm/plot.index,data=concatenate,method='ML')

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

#best intercepts model:
modC1<-lme(sqrt(spalpha.rich)~Ca,~1|comm/plot.index,data=concatenate,method='ML')
diagnostics(modC1)
#df: 85
#AICc: 202.5
#marginal-r2: 9%
#conditional-r2: 88%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(sqrt(spalpha.rich)~1,~predictor-1|comm/plot.index,control=ctrl,data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

#best slopes model:
modC2<-lme(sqrt(spalpha.rich)~1,~Cappm-1|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modC2)
#df: 85
#AICc: 272.9
#marginal-r2: 0%
#conditional-r2: 86%

#D1
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(sqrt(spalpha.rich)~Ca,~predictor|comm/plot.index,control=ctrl,data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

#best model constrained by intercepts:
modD1<-lme(sqrt(spalpha.rich)~Ca,~spalpha.light|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modD1)
#df: 1, 84
#AICc: 206.1
#marginal-r2: 8%
#conditional-r2: 89%
modD1.exp<-lme(sqrt(spalpha.rich)~Ca+spalpha.light,~spalpha.light|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 2, 83
#AICc: 208.3
#marginal-r2: 9%
#conditional-r2: 89%

#D2
global.model<-lme(sqrt(spalpha.rich) ~ OM + P + Kppm + Mgppm + Cappm + 
                    pH + CEC + K + Mg + Ca + C + N + C.N + 
                    spalpha.light + spalpha.l + spalpha.o + 
                    spalpha.r, ~Cappm|comm/plot.index,control=ctrl,data=concatenate,method='ML')
model<-dredge(global.model,subset=smat)#,subset = smat,extra = list(
  #"R^2", "*" = function(x) {
   # c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  #}))
models<-model[model$delta<2,]
plot(models)

#best model constrained by slopes:
modD2<-lme(sqrt(spalpha.rich)~Ca,~Cappm|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modD2)
#df: 1, 84
#AICc: 209.9
#marginal-r2: 11%
#conditional-r2: 87%
modD2.exp<-lme(sqrt(spalpha.rich)~Ca+Cappm,~Cappm|comm/plot.index,control=ctrl,data=concatenate,method='ML')
diagnostics(modD2.exp)
#df: 2, 83
#AICc: 210.9
#marginal-r2: 13%
#conditional-r2: 88%
########################################################################
#plot best model(s):
########################################################################
#plot best subplots in plots model modC1 
plot(sqrt(spalpha.rich)~Ca,data=concatenate,pch=NA,
     main='plot model C1',
     ylab='subplot-level richness',xlab='normalized calcium',axes=FALSE)
text(concatenate$Ca,sqrt(concatenate$spalpha.rich),labels=concatenate$plot.index,cex=1)
axis(1,at=c(-2,-1.5,-1,-.5,0,.5,1,1.5,2),labels=c(-2,-1.5,-1,-.5,0,.5,1,1.5,2),cex.axis=1.5)
axis(2,at=c(0,1,2,3,4,5),labels=round((c(0,1,2,3,4,5))^2,1),cex.axis=1.5)
box()
abline(2.6168275,0.3734736,col='red',lwd=5)

for(i in 1:dim(coef(modC1))[1]){abline(coef(modC1)[i,1],0.3734736,lwd=1)}
legend('topleft',bty='n',legend=c('across subplots','within plots'),col=c('red','black'),lwd=c(5,1),pch=c(NA,NA),cex=1.5)

########################################################################
#ALPHA EVENNESS
#(4)
########################################################################
concatenate<-data.frame(spalpha.ab,spalpha.patterns.f)
concatenate<-concatenate[c(1:90,95:112,114),]
rich<-scale(concatenate$spalpha.rich,center=TRUE)
#data<-data.frame(rich,spalpha.ab[c(1:90,95:112,114),2:18])

data<-spalpha.ab[c(1:90,95:112,114),c(2:4,8:10,12:18)]
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .4, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
smat <- outer(1:dim(data)[2], 1:dim(data)[2], vCorrelated, data = data)
nm <- colnames(data)
dimnames(smat) <- list(nm, nm)
#don't care about Mgppm, Cappm, pH, Ca, because they're correlated with richness

########################################################################
#modeling subplots in plots
########################################################################
#step (B)
modB<-lme(log(spalpha.evens)~1,random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modB)
#df: 81
#AICc: 201.8
#marginal-r2: 0%
#conditional-r2: 51%

#step (C1)
global.model<-lme(log(spalpha.evens) ~ OM + P + Kppm + CEC + K + Mg + C +N+C.N+ 
                    spalpha.light + spalpha.l + spalpha.o + 
                    spalpha.r, random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

#best intercepts model
modC1<-lme(log(spalpha.evens)~Mg,random=list(plot.index=pdDiag(~1)),data=concatenate,method='ML')
diagnostics(modC1)
#df: 1, 80
#AICc: 201.1
#marginal-r2: 3%
#conditional-r2: 49%

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(log(spalpha.evens)~1,random=list(plot.index=pdDiag(~predictor-1)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

#best slopes model:
modC2<-lme(log(spalpha.evens)~1,random=list(plot.index=pdDiag(~OM-1)),data=concatenate,method='ML')
diagnostics(modC2)
#df: 81
#AICc: 218.0
#marginal-r2: 0%
#conditional-r2: 36%

#step (D1)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(log(spalpha.evens)~Mg,random=list(plot.index=pdDiag(~predictor)),data=concatenate,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

#best model constrained by intercepts
modD1<-lme(log(spalpha.evens)~Mg,random=list(plot.index=pdDiag(~C.N)),data=concatenate,method='ML')
diagnostics(modD1)
#df: 1, 80
#AICc: 202.2
#marginal-r2: 1%
#conditional-r2: 60%
modD1.exp<-lme(log(spalpha.evens)~Mg+C.N,random=list(plot.index=pdDiag(~C.N)),data=concatenate,method='ML')
diagnostics(modD1.exp)
#df: 2, 79
#AICc: 204.2
#marginal-r2: 2%
#conditional-r2: 58%

#step (D2)
global.model<-lme(log(spalpha.evens) ~ OM + P + Kppm + CEC + K + Mg + C +N+C.N+ 
                    spalpha.light + spalpha.l + spalpha.o + 
                    spalpha.r, random=list(plot.index=pdDiag(~OM)),data=concatenate,method='ML')

model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

#best model constrained by slopes
modD2<-lme(log(spalpha.evens)~Mg,random=list(plot.index=pdDiag(~OM)),data=concatenate,method='ML')
diagnostics(modD2)
#df: 1, 80
#AICc: 202.4
#marginal-r2: 3%
#conditional-r2: 53%
modD2.exp<-lme(log(spalpha.evens)~Mg+OM,random=list(plot.index=pdDiag(~OM)),data=concatenate,method='ML')
diagnostics(modD2.exp)
#df: 2, 79
#AICc: 204.6
#marginal-r2: 3%
#conditional-r2: 53%
########################################################################
#modeling subplots in plots nested in forests
########################################################################
#step (B)
modB<-lme(log(spalpha.evens)~1,~1|comm/plot.index,data=concatenate,method='ML')
diagnostics(modB)
#df: 81
#AICc: 201.8
#marginal-r2: 0%
#conditional-r2: 51%

#step (C1)
global.model<-lme(log(spalpha.evens) ~ OM + P + Kppm + CEC + K + Mg + C + N+C.N+
                    spalpha.light + spalpha.l + spalpha.o + 
                    spalpha.r, ~1|comm/plot.index,data=concatenate,method='ML')
model<-dredge(global.model,subset = smat,extra = list(
  "R^2", "*" = function(x) {
    c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  }))
models<-model[model$delta<2,]
plot(models)

modC1<-lme(log(spalpha.evens)~1,~1|comm/plot.index,data=concatenate,method='ML')
diagnostics(modC1)
#same as unconditional model

#step (C2)
score<-NULL
for(i in 1:length(colnames(concatenate[2:18]))){
  predictor<-concatenate[,i+1]
  model<-lme(log(spalpha.evens)~1,~predictor-1|comm/plot.index,data=concatenate,control=ctrl,method='ML')
  score[i]<-AICc(model)
}
results<-data.frame(colnames(concatenate[,2:18]),score)
View(results[order(results[,2]),])

#best slopes model
modC2<-lme(log(spalpha.evens)~1,~spalpha.r-1|comm/plot.index,data=concatenate,control=ctrl,method='ML')
diagnostics(modC2)
#df: 81
#AICc: 219.1
#marginal-r2: 0%
#conditional-r2: 31%

#step (D1)
#by hand:
model<-lme(log(spalpha.evens)~1,~predictor|comm/plot.index,data=concatenate,control=ctrl,method='ML')
AICc(model)

#best model constrained by intercepts
modD1<-lme(log(spalpha.evens)~1,~spalpha.r|comm/plot.index,data=concatenate,control=ctrl,method='ML')
diagnostics(modD1)
#df: 81
#AICc: 206.4
#marginal-r2: 0%
#conditional-r2: 52%
modD1.exp<-lme(log(spalpha.evens)~spalpha.r,~spalpha.r|comm/plot.index,data=concatenate,control=ctrl,method='ML')
diagnostics(modD1.exp)
#df: 81
#AICc: 208.5
#marginal-r2: 0%
#conditional-r2: 52%

#step (D2)
global.model<-lme(log(spalpha.evens) ~ OM + P + Kppm + CEC + K + Mg + C + N+C.N+
                    spalpha.light + spalpha.l + spalpha.o + 
                    spalpha.r, ~spalpha.r|comm/plot.index,control=ctrl,data=concatenate,method='ML')
model<-dredge(global.model,subset = smat)#,extra = list(
  #"R^2", "*" = function(x) {
   # c(MRsq = r.squaredGLMM(x)[1], CRsq = r.squaredGLMM(x)[2])
  #}))
models<-model[model$delta<2,]
plot(models)

#best model constrained by slopes:
modD2<-lme(log(spalpha.evens)~1,~spalpha.r|comm/plot.index,data=concatenate,control=ctrl,method='ML')
diagnostics(modD2)
#same as D1
########################################################################
#plot best model(s)
########################################################################
#plot model with the lowest AICc score modC1 with subplots in plots
plot(log(spalpha.evens)~Mg,data=concatenate,pch=NA,
     main='Spring Forbs: abiotic conditions 
     and subplot-level evenness',ylab='subplot-level evenness',
     xlab='normalized Mg',yaxt="n",cex.axis=1.5,xlim=c(-1.75,3.75))
#xlim=c(-1.75,4))
text(concatenate$Mg,log(concatenate$spalpha.evens),labels=concatenate$plot.index,cex=1)
axis(2,at=c(-2,-1.5,-1,-.5,0),las=2,labels=round(c(exp(-2),exp(-1.5),exp(-1),exp(-.5),exp(0)),digits=1),cex.axis=1.5)
for(i in 1:dim(coef(modC1))[1]){
  #oldline<-(coef(modC1)[i,1]+((coef(modC1)[i,2])*(concatenate$C.N[concatenate$plot.index==i])))
  #newslope<-coef(modC1)[i,2]*coef(lm(concatenate$spalpha.evens[concatenate$plot.index==i]~oldline))[2]
  abline(coef(modC1)[i,1],-0.1280336)}
abline(-0.9431351,-0.1280336,col='red',lwd=5)
legend('bottomright',cex=1.5,legend=c('across subplots','within plots'),col=c('red','black'),lwd=c(3,1))

#plot model with the highest variance explained modD1 with subplots in plots:
plot(log(spalpha.evens)~Mg,data=concatenate,pch=NA,
     main='Spring Forbs: abiotic conditions 
     and subplot-level evenness',ylab='subplot-level evenness',
     xlab='normalized Mg',yaxt="n")
#xlim=c(-1.75,4))
text(concatenate$Mg,log(concatenate$spalpha.evens),labels=concatenate$plot.index,cex=.7)
axis(2,at=c(-2,-1.5,-1,-.5,0),las=2,labels=round(c(exp(-2),exp(-1.5),exp(-1),exp(-.5),exp(0)),digits=1))
for(i in 1:dim(coef(modD1))[1]){abline(coef(modD1)[i,1],-0.09327729)}
abline(-0.9472909,-0.0932773,col='red',lwd=3)
legend('bottomright',legend=c('across subplots','within plots'),col=c('red','black'),lwd=c(3,1))

plot(log(spalpha.evens)~C.N,data=concatenate,pch=NA,
     main='Spring Forbs: abiotic conditions 
     and subplot-level evenness',ylab='subplot-level evenness',
     xlab='normalized C:N',yaxt="n",cex.axis=1.5)
#,ylim=c(0,1))
#xlim=c(-1.75,4))
text(concatenate$C.N,log(concatenate$spalpha.evens),labels=concatenate$plot.index,cex=1)
axis(2,at=c(-2,-1.5,-1,-.5,0),las=2,labels=round(c(exp(-2),exp(-1.5),exp(-1),exp(-.5),exp(0)),digits=1),cex.axis=1.5)
for(i in 1:dim(coef(modD1))[1]){
  #oldline<-(coef(modD1)[i,1]+((coef(modD1)[i,2])*(concatenate$C.N[concatenate$plot.index==i])))
  #newslope<-coef(modD1)[i,2]*coef(lm(concatenate$spalpha.evens[concatenate$plot.index==i]~oldline))[2]
  abline(coef(modD1)[i,1],coef(modD1)[i,3])}
