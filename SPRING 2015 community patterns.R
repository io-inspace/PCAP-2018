#SPRING 2015 COMMUNITY PATTERNS
setwd("C:/Users/Sam/Desktop/PCAP")

#call in modified (from the unweildy PCAP sheets) spring data, and the 
#seasonal designation data
springbymod<-read.csv("SPRING 2015 response.csv",header=TRUE,sep=',')
desig<-read.csv("SPRING 2015 seasonal designations.csv",header=TRUE,sep=',')

designations<-1 #1 to remove summer forbs

#if you want to add plots and mods:
#for subplots:
plot<-as.vector(c(rep(c(1:5),each=4),rep(6,2),rep(7:29,each=4)))
mods<-as.vector(c(rep(c('a','b','c','d'),5),'a','b',rep(c('a','b','c','d'),23)))
plots.mods<-data.frame(plot,mods)
#for plots:
plots<-c(1:29)

#get a list of forms in spring by mod
for(i in 1:dim(springbymod)[1]){
  for(j in 1:dim(desig)[1]){
    if(springbymod$genus[i]==desig$genus[j] && 
       springbymod$species[i]==desig$species[j])
    {springbymod$form[i]<-as.character(desig$form[j])}
  }
}

#get seasonal designations
for(i in 1:dim(springbymod)[1]){
  for(j in 1:dim(desig)[1]){
    if(springbymod$genus[i]==desig$genus[j] && 
       springbymod$species[i]==desig$species[j])
    {springbymod$month[i]<-as.character(desig$designation[j])}
  }
}

#REMOVE SUMMER FORBS?
if(designations == 1){
springspecies<-springbymod[!springbymod$month == "summer",]}else{
  springspecies<-springbymod
}

#REMOVE NON-FORBS?
#springspecies<-springspecies[springspecies$form == "forb",]

#RUN FROM HERE FOR SUBSEQUENT VERSIONS

otherforms<-0 #to remove forbs
forbsonly<-1 #1 to remove non-forbs
z<-0 #to return scaled values
normalize<-0 #1 to return normalized (log- or sqrt-transformed data)
corrstructure<-0 #1 to return visualization

if(otherforms == 1){
  
springspecies<-springspecies[springspecies$form != "forb",]}
  
  
#richness

#alpha
spmods.rich<-split(springspecies[,c(7,11,15,19)],springspecies[,1],drop=FALSE)
spalpha.rich<-matrix(data=NA,ncol=dim(spmods.rich[[1]])[2],nrow=length(spmods.rich))
for(i in 1:length(spmods.rich)){for(j in 1:dim(spmods.rich[[i]])[2]){
  spalpha.rich[i,j]<-sum(spmods.rich[[i]][j])
}
}
spalpha.rich[spalpha.rich == 0]<-NA
colnames(spalpha.rich)<-c("subplot1.rich","subplot2.rich","subplot3.rich",
                          "subplot4.rich")

spalpha.rich.keep<-spalpha.rich #keep a version in this format


#mean alpha
spmean.rich<-apply(spalpha.rich, 1, mean,na.rm=TRUE)

#gamma
spgamma.rich<-NULL
for(i in 1:length(spmods.rich)){spgamma.rich[i]<-dim(spmods.rich[[i]])[1]}

#beta
spbeta.rich<-(spmean.rich/spgamma.rich)

#variance
library(raster)
spvar.rich<-apply(spalpha.rich,1,cv,na.rm=TRUE)

#richness dataframes
#across subplots:

spalpha.rich<-as.matrix(unlist(as.list(t(spalpha.rich))))
spalpha.rich<-spalpha.rich[!is.na(spalpha.rich)]

#across plots:
spplot.rich<-data.frame(spgamma.rich,spmean.rich,
                        spbeta.rich,spvar.rich)

#NORMALIZE?
if(normalize==1){
  spgamma.rich.n<-sqrt(spgamma.rich)
  spmean.rich.n<-sqrt(spmean.rich)
  spbeta.rich.n<-sqrt(spbeta.rich)
  spvar.rich.n<-log(spvar.rich)
  spplot.rich.n<-data.frame(spgamma.rich.n,spmean.rich.n,
                            spbeta.rich.n,spvar.rich.n)
  spalpha.rich.n<-sqrt(spalpha.rich)
}


#evenness
midpoints<-c(.025,.05,1.5,3.5,7.5,17.5,37.5,62.5,85,97.5)
spmods.cover<-split(springspecies[,c(6,10,14,18)],springspecies[,1],drop=FALSE)

for(i in 1:length(spmods.cover)){
  for(j in 1:dim(spmods.cover[[i]])[2]){
    for(k in 1:dim(spmods.cover[[i]])[1]){if(is.na(spmods.cover[[i]][k,j])){spmods.cover[[i]][k,j]<-0}
      else(spmods.cover[[i]][k,j]<-midpoints[spmods.cover[[i]][k,j]])}
  }
}
spmods.cover.keep<-spmods.cover

spmods.prop<-spmods.cover
for(i in 1:length(spmods.cover)){
  for(j in 1:dim(spmods.cover[[i]])[2]){
    for(k in 1:dim(spmods.cover[[i]])[1]){
      total<-sum(spmods.cover[[i]][,j])
      spmods.prop[[i]][k,j]<-(spmods.cover[[i]][k,j]/total)^2}
  }
}

spalpha.evens<-matrix(data=NA,ncol=dim(spmods.cover[[1]])[2],
                      nrow=length(spmods.cover),byrow=TRUE)
for(i in 1:length(spmods.cover)){for(j in 1:dim(spmods.cover[[i]])[2]){
  spalpha.evens[i,j]<-(1/sum(spmods.prop[[i]][,j]))
}
}

colnames(spalpha.evens)<-c("subplot1.evens","subplot2.evens","subplot3.evens",
                           "subplot4.evens")

spalpha.evens<-(spalpha.evens/spalpha.rich.keep)

spalpha.evens.keep<-spalpha.evens #keep a version in this format

#mean alpha
spmean.evens<-apply(spalpha.evens, 1, mean,na.rm=TRUE)

#gamma
spplots.cover<-list()
for(i in 1:length(spmods.cover.keep)){
  spplots.cover[[i]]<-apply(spmods.cover.keep[[i]],1,sum)
}

spgamma.evens<-NULL
for(i in 1:length(spplots.cover)){
  total<-sum(spplots.cover[[i]])
  proportions<-NULL
  for(j in 1:length(spplots.cover[[i]])){
    proportions[j]<-(spplots.cover[[i]][j]/total)^2}
  spgamma.evens[i]<-1/sum(proportions)
}

spgamma.evens<-(spgamma.evens/spgamma.rich)

#beta
spbeta.evens<-(spmean.evens/spgamma.evens) 
#remember beta here is not on a 0-1 scale; if it's less than 1, that means
#that gamma evenness > mean evenness.  If it's greater than 1, that means that
#gamma evenness < mean evenness. Interpret carefully because it's weird.

#variance
spvar.evens<-apply(spalpha.evens,1,cv,na.rm=TRUE)

#evenness dataframes
#plot
spplot.evens<-data.frame(spgamma.evens,spmean.evens,
                         spbeta.evens,spvar.evens)

#subplot
spalpha.evens<-as.matrix(unlist(as.list(t(spalpha.evens))))
spalpha.evens<-spalpha.evens[!is.na(spalpha.evens)]

if(normalize==1){
  spgamma.evens.n<-log(spgamma.evens)
  spmean.evens.n<-log(spmean.evens)
  spbeta.evens.n<-log(spbeta.evens)
  spvar.evens.n<-sqrt(spvar.evens)
  spplot.evens.n<-data.frame(spgamma.evens.n, spmean.evens.n,
                             spbeta.evens.n,spvar.evens.n)
  spalpha.evens.n<-log(spalpha.evens)
}

#all patterns
spplot.patterns<-data.frame(spplot.rich,spplot.evens)
spalpha.patterns<-data.frame(spalpha.rich,spalpha.evens)

#NORMALIZE?
if(normalize==1){
  spplot.patterns.n<-data.frame(spplot.rich.n,
                                spplot.evens.n)
  spalpha.patterns.n<-data.frame(spalpha.rich.n,spalpha.evens.n)
  if(z == 1){
    spplot.patterns.n.z<-scale(spplot.patterns.n, center=TRUE,scale=TRUE)
    spalpha.patterns.n.z<-scale(spalpha.patterns.n,center=TRUE,scale=TRUE)
  }
}

if(z == 1){
  spplot.patterns.z<-scale(spplot.patterns, center=TRUE,scale=TRUE)
  spalpha.patterns.z<-scale(spalpha.patterns,center=TRUE,scale=TRUE)
}

#save versions:
if(z==0 & normalize == 0){
  spplot.patterns.otherforms<-spplot.patterns
  spalpha.patterns.otherforms<-spalpha.patterns
}else
  if(z==0 & normalize == 1){
    spplot.patterns.n.otherforms<-spplot.patterns.n
    spalpha.patterns.n.otherforms<-spalpha.patterns.n
  }else
    if(z==1 & normalize == 0){
      spplot.patterns.z.otherforms<-spplot.patterns.z
      spalpha.patterns.z.otherforms<-spalpha.patterns.z
    }else
      if(z==1 & normalize==1){
        spplot.patterns.n.z.otherforms<-spplot.patterns.n.z
        spalpha.patterns.n.z.otherforms<-spalpha.patterns.n.z
      }

#VISUALIZE?
if(corrstructure==1){
  library(PerformanceAnalytics)
  chart.Correlation(spplot.rich.n)
  chart.Correlation(spplot.evens.n)
  chart.Correlation(spplot.patterns.n)
  chart.Correlation(spalpha.patterns.n)
}


if(forbsonly == 0){

#richness

#alpha
spmods.rich<-split(springspecies[,c(7,11,15,19)],springspecies[,1],drop=FALSE)
spalpha.rich<-matrix(data=NA,ncol=dim(spmods.rich[[1]])[2],nrow=length(spmods.rich))
for(i in 1:length(spmods.rich)){for(j in 1:dim(spmods.rich[[i]])[2]){
  spalpha.rich[i,j]<-sum(spmods.rich[[i]][j])
}
}
spalpha.rich[spalpha.rich == 0]<-NA
colnames(spalpha.rich)<-c("subplot1.rich","subplot2.rich","subplot3.rich",
                          "subplot4.rich")

spalpha.rich.keep<-spalpha.rich #keep a version in this format


#mean alpha
spmean.rich<-apply(spalpha.rich, 1, mean,na.rm=TRUE)

#gamma
spgamma.rich<-NULL
for(i in 1:length(spmods.rich)){spgamma.rich[i]<-dim(spmods.rich[[i]])[1]}

#beta
spbeta.rich<-(spmean.rich/spgamma.rich)

#variance
library(raster)
spvar.rich<-apply(spalpha.rich,1,cv,na.rm=TRUE)

#richness dataframes
#across subplots:

spalpha.rich<-as.matrix(unlist(as.list(t(spalpha.rich))))
spalpha.rich<-spalpha.rich[!is.na(spalpha.rich)]

#across plots:
spplot.rich<-data.frame(spgamma.rich,spmean.rich,
                        spbeta.rich,spvar.rich)

#NORMALIZE?
if(normalize==1){
  spgamma.rich.n<-sqrt(spgamma.rich)
  spmean.rich.n<-sqrt(spmean.rich)
  spbeta.rich.n<-sqrt(spbeta.rich)
  spvar.rich.n<-log(spvar.rich)
  spplot.rich.n<-data.frame(spgamma.rich.n,spmean.rich.n,
                               spbeta.rich.n,spvar.rich.n)
  spalpha.rich.n<-sqrt(spalpha.rich)
}


#evenness
midpoints<-c(.025,.05,1.5,3.5,7.5,17.5,37.5,62.5,85,97.5)
spmods.cover<-split(springspecies[,c(6,10,14,18)],springspecies[,1],drop=FALSE)

for(i in 1:length(spmods.cover)){
  for(j in 1:dim(spmods.cover[[i]])[2]){
    for(k in 1:dim(spmods.cover[[i]])[1]){if(is.na(spmods.cover[[i]][k,j])){spmods.cover[[i]][k,j]<-0}
      else(spmods.cover[[i]][k,j]<-midpoints[spmods.cover[[i]][k,j]])}
  }
}
spmods.cover.keep<-spmods.cover

spmods.prop<-spmods.cover
for(i in 1:length(spmods.cover)){
  for(j in 1:dim(spmods.cover[[i]])[2]){
    for(k in 1:dim(spmods.cover[[i]])[1]){
      total<-sum(spmods.cover[[i]][,j])
      spmods.prop[[i]][k,j]<-(spmods.cover[[i]][k,j]/total)^2}
  }
}

spalpha.evens<-matrix(data=NA,ncol=dim(spmods.cover[[1]])[2],
                      nrow=length(spmods.cover),byrow=TRUE)
for(i in 1:length(spmods.cover)){for(j in 1:dim(spmods.cover[[i]])[2]){
  spalpha.evens[i,j]<-(1/sum(spmods.prop[[i]][,j]))
  }
}

colnames(spalpha.evens)<-c("subplot1.evens","subplot2.evens","subplot3.evens",
                        "subplot4.evens")

spalpha.evens<-(spalpha.evens/spalpha.rich.keep)

spalpha.evens.keep<-spalpha.evens #keep a version in this format

#mean alpha
spmean.evens<-apply(spalpha.evens, 1, mean,na.rm=TRUE)

#gamma
spplots.cover<-list()
for(i in 1:length(spmods.cover.keep)){
  spplots.cover[[i]]<-apply(spmods.cover.keep[[i]],1,sum)
}

spgamma.evens<-NULL
for(i in 1:length(spplots.cover)){
  total<-sum(spplots.cover[[i]])
  proportions<-NULL
  for(j in 1:length(spplots.cover[[i]])){
    proportions[j]<-(spplots.cover[[i]][j]/total)^2}
  spgamma.evens[i]<-1/sum(proportions)
}

spgamma.evens<-(spgamma.evens/spgamma.rich)

#beta
spbeta.evens<-(spmean.evens/spgamma.evens) 
#remember beta here is not on a 0-1 scale; if it's less than 1, that means
#that gamma evenness > mean evenness.  If it's greater than 1, that means that
#gamma evenness < mean evenness. Interpret carefully because it's weird.
                                    
#variance
spvar.evens<-apply(spalpha.evens,1,cv,na.rm=TRUE)

#evenness dataframes
#plot
spplot.evens<-data.frame(spgamma.evens,spmean.evens,
                         spbeta.evens,spvar.evens)

#subplot
spalpha.evens<-as.matrix(unlist(as.list(t(spalpha.evens))))
spalpha.evens<-spalpha.evens[!is.na(spalpha.evens)]

if(normalize==1){
  spgamma.evens.n<-log(spgamma.evens)
  spmean.evens.n<-log(spmean.evens)
  spbeta.evens.n<-log(spbeta.evens)
  spvar.evens.n<-sqrt(spvar.evens)
  spplot.evens.n<-data.frame(spgamma.evens.n, spmean.evens.n,
                                spbeta.evens.n,spvar.evens.n)
  spalpha.evens.n<-log(spalpha.evens)
}

#all patterns
spplot.patterns<-data.frame(spplot.rich,spplot.evens)
spalpha.patterns<-data.frame(spalpha.rich,spalpha.evens)

#NORMALIZE?
if(normalize==1){
  spplot.patterns.n<-data.frame(spplot.rich.n,
                                   spplot.evens.n)
  spalpha.patterns.n<-data.frame(spalpha.rich.n,spalpha.evens.n)
  if(z == 1){
    spplot.patterns.n.z<-scale(spplot.patterns.n, center=TRUE,scale=TRUE)
    spalpha.patterns.n.z<-scale(spalpha.patterns.n,center=TRUE,scale=TRUE)
  }
}

if(z == 1){
  spplot.patterns.z<-scale(spplot.patterns, center=TRUE,scale=TRUE)
  spalpha.patterns.z<-scale(spalpha.patterns,center=TRUE,scale=TRUE)
}

#save versions:
if(z==0 & normalize == 0){
  spplot.patterns.allforms<-spplot.patterns
  spalpha.patterns.allforms<-spalpha.patterns
}else
if(z==0 & normalize == 1){
  spplot.patterns.n.allforms<-spplot.patterns.n
  spalpha.patterns.n.allforms<-spalpha.patterns.n
}else
if(z==1 & normalize == 0){
  spplot.patterns.z.allforms<-spplot.patterns.z
  spalpha.patterns.z.allforms<-spalpha.patterns.z
}else
if(z==1 & normalize==1){
  spplot.patterns.n.z.allforms<-spplot.patterns.n.z
  spalpha.patterns.n.z.allforms<-spalpha.patterns.n.z
}

#VISUALIZE?
if(corrstructure==1){
  library(PerformanceAnalytics)
  chart.Correlation(spplot.rich.n)
  chart.Correlation(spplot.evens.n)
  chart.Correlation(spplot.patterns.n)
  chart.Correlation(spalpha.patterns.n)
  }
}

############################################################################

#REMOVE NON-FORBS?
if(forbsonly == 1){springspecies<-springspecies[springspecies$form == "forb",]

#richness

#alpha
spmods.rich<-split(springspecies[,c(7,11,15,19)],springspecies[,1],drop=FALSE)
spalpha.rich<-matrix(data=NA,ncol=dim(spmods.rich[[1]])[2],nrow=length(spmods.rich))
for(i in 1:length(spmods.rich)){for(j in 1:dim(spmods.rich[[i]])[2]){
  spalpha.rich[i,j]<-sum(spmods.rich[[i]][j])
}
}
spalpha.rich[spalpha.rich == 0]<-NA
colnames(spalpha.rich)<-c("subplot1.rich","subplot2.rich","subplot3.rich",
                          "subplot4.rich")
#there are no forbs in 24.
twentyfour<-matrix(rep(0,4),ncol=4,nrow=1)
colnames(twentyfour)<-c("subplot1.rich","subplot2.rich","subplot3.rich",
                        "subplot4.rich")
spalpha.rich<-rbind(spalpha.rich[1:23,],twentyfour,spalpha.rich[24:28,])

spalpha.rich.keep<-spalpha.rich #keep a version in this format

#mean alpha
spmean.rich<-apply(spalpha.rich, 1, mean,na.rm=TRUE)

#gamma
spgamma.rich<-NULL
for(i in 1:length(spmods.rich)){spgamma.rich[i]<-dim(spmods.rich[[i]])[1]}
spgamma.rich<-c(spgamma.rich[1:23],0,spgamma.rich[24:28])

#beta
spbeta.rich<- (spmean.rich/spgamma.rich)
spbeta.rich[spbeta.rich=="NaN"]<-NA

#variance
library(raster)
spvar.rich<-apply(spalpha.rich,1,cv,na.rm=TRUE)

#richness dataframes
#across subplots:
spalpha.rich<-as.matrix(unlist(as.list(t(spalpha.rich))))
spalpha.rich[115,]<-0
spalpha.rich<-spalpha.rich[!is.na(spalpha.rich)]

#across plots:
spplot.rich<-data.frame(spgamma.rich,spmean.rich,
                        spbeta.rich,spvar.rich)

#NORMALIZE?
if(normalize==1){
  spmean.rich[spvar.rich == 0]<-.000001
  spmean.rich.n<-sqrt(spmean.rich)
  spgamma.rich[spgamma.rich == 0]<-.00001
  spgamma.rich.n<-sqrt(spgamma.rich)
  spvar.rich[spvar.rich == 0]<-.000001
  spvar.rich.n<-sqrt(spvar.rich)
  spplot.rich.n<-data.frame(spgamma.rich.n,spmean.rich.n,
                               spbeta.rich,spvar.rich.n)
  spalpha.rich.n<-sqrt(spalpha.rich)
}

#evenness
midpoints<-c(.025,.05,1.5,3.5,7.5,17.5,37.5,62.5,85,97.5)
spmods.cover<-split(springspecies[,c(6,10,14,18)],springspecies[,1],drop=FALSE)

for(i in 1:length(spmods.cover)){
  for(j in 1:dim(spmods.cover[[i]])[2]){
    for(k in 1:dim(spmods.cover[[i]])[1]){if(is.na(spmods.cover[[i]][k,j])){spmods.cover[[i]][k,j]<-0}
      else(spmods.cover[[i]][k,j]<-midpoints[spmods.cover[[i]][k,j]])}
  }
}
spmods.cover.keep<-spmods.cover

spmods.prop<-spmods.cover
for(i in 1:length(spmods.cover)){
  for(j in 1:dim(spmods.cover[[i]])[2]){
    for(k in 1:dim(spmods.cover[[i]])[1]){
      total<-sum(spmods.cover[[i]][,j])
      spmods.prop[[i]][k,j]<-(spmods.cover[[i]][k,j]/total)^2}
  }
}

spalpha.evens<-matrix(data=NA,ncol=dim(spmods.cover[[1]])[2],
                      nrow=length(spmods.cover),byrow=TRUE)
for(i in 1:length(spmods.cover)){for(j in 1:dim(spmods.cover[[i]])[2]){
  spalpha.evens[i,j]<-(1/sum(spmods.prop[[i]][,j]))
}
}

twentyfour<-matrix(rep(NA,4),ncol=4,nrow=1)
colnames(twentyfour)<-c("subplot1.rich","subplot2.rich","subplot3.rich",
                        "subplot4.rich")
spalpha.evens<-rbind(spalpha.evens[1:23,],twentyfour,spalpha.evens[24:28,])

spalpha.evens.keep<-spalpha.evens #keep a version in this format

spalpha.evens<-(spalpha.evens/spalpha.rich.keep)

#mean alpha
spmean.evens<-apply(spalpha.evens, 1, mean,na.rm=TRUE)
spmean.evens[spmean.evens=="NaN"]<-NA

#gamma
spplots.cover<-list()
for(i in 1:length(spmods.cover.keep)){
  spplots.cover[[i]]<-apply(spmods.cover.keep[[i]],1,sum)
}

spgamma.evens<-NULL
for(i in 1:length(spplots.cover)){
  total<-sum(spplots.cover[[i]])
  proportions<-NULL
  for(j in 1:length(spplots.cover[[i]])){
    proportions[j]<-(spplots.cover[[i]][j]/total)^2}
  spgamma.evens[i]<-1/sum(proportions)
}
spgamma.evens<-c(spgamma.evens[1:23],NA,spgamma.evens[24:28])

spgamma.evens<-(spgamma.evens/spgamma.rich)

#beta
spbeta.evens<-(spmean.evens/spgamma.evens)

#variance
spvar.evens<-apply(spalpha.evens,1,cv,na.rm=TRUE)

#evenness dataframes
#plot
spplot.evens<-data.frame(spgamma.evens,spmean.evens,
                         spbeta.evens,spvar.evens)

#subplot
spalpha.evens<-as.matrix(unlist(as.list(t(spalpha.evens))))
spalpha.evens<-c(spalpha.evens[1:22],spalpha.evens[25:116])
spalpha.evens[spalpha.evens=="NaN"]<-NA

if(normalize==1){
  spgamma.evens[spgamma.evens == 0]<-.000001
  spgamma.evens.n<-log(spgamma.evens)
  spmean.evens[spmean.evens == 0]<-.000001
  spmean.evens.n<-sqrt(spmean.evens)
  spbeta.evens[spbeta.evens == 0]<-.000001
  spbeta.evens.n<-log(spbeta.evens) #not really normalized, too much left skew.
  spvar.evens[spvar.evens == 0]<-.000001
  spvar.evens.n<-sqrt(spvar.evens)
  spplot.evens.n<-data.frame(spgamma.evens.n, spmean.evens.n,
                                spbeta.evens.n,spvar.evens.n)
  spalpha.evens.n<-log(spalpha.evens)
}

#all patterns
spplot.patterns<-data.frame(spplot.rich,spplot.evens)
spalpha.patterns<-data.frame(spalpha.rich,spalpha.evens)
if(z == 1){
  spplot.patterns.z<-scale(spplot.patterns,center=TRUE,scale=TRUE)
  spalpha.patterns.z<-scale(spalpha.patterns,center=TRUE,scale=TRUE)
}
#NORMALIZE?
if(normalize==1){
  spplot.patterns.n<-data.frame(spplot.rich.n,
                                   spplot.evens.n)
  spalpha.patterns.n<-data.frame(spalpha.rich.n,spalpha.evens.n)
  if(z == 1){
    spplot.patterns.n.z<-scale(spplot.patterns.n,center=TRUE,scale=TRUE)
    spalpha.patterns.n.z<-scale(spalpha.patterns.n,center=TRUE,scale=TRUE)
  }
}

#save versions:
if(z==0 & normalize == 0){
  spplot.patterns.forbs<-spplot.patterns
  spalpha.patterns.forbs<-spalpha.patterns
}else
  if(z==0 & normalize == 1){
    spplot.patterns.n.forbs<-spplot.patterns.n
    spalpha.patterns.n.forbs<-spalpha.patterns.n
  }else
    if(z==1 & normalize == 0){
      spplot.patterns.z.forbs<-spplot.patterns.z
      spalpha.patterns.z.forbs<-spalpha.patterns.z
    }else
      if(z==1 & normalize==1){
        spplot.patterns.n.z.forbs<-spplot.patterns.n.z
        spalpha.patterns.n.z.forbs<-spalpha.patterns.n.z
      }

#VISUALIZE?
if(corrstructure==1){
  library(PerformanceAnalytics)
  chart.Correlation(spplot.rich.n)
  chart.Correlation(spplot.evens.n)
  chart.Correlation(spplot.patterns.n)
  chart.Correlation(spalpha.patterns.n)
  }
}
