#SUMMER 2015 COMMUNITY PATTERNS
setwd("C:/Users/Sam/Desktop/PCAP")

#call in modified (from the unweildy PCAP sheets) summer data
summerbymod<-read.csv("SUMMER 2015 response.csv",header=TRUE,sep=',')
desig<-read.csv("SUMMER 2015 seasonal designations_edited.csv",header=TRUE,sep=',')

#remove any species not present in the herbaceous layer:
summerbymod<-summerbymod[summerbymod$herb.y.n == 1,]
#remove species only present in residuals:
summerbymod<-summerbymod[summerbymod$res == 0,]

#makes some tough choices:
#binary 1 for yes, 0 for no
designations<-0 #1 to remove spring forbs

#if you want to add plots and mods:
#for subplots:
plot<-as.vector(c(rep(c(1:5),each=4),rep(6,2),rep(7:29,each=4)))
mods<-as.vector(c(rep(c('a','b','c','d'),5),'a','b',rep(c('a','b','c','d'),23)))
plots.mods<-data.frame(plot,mods)
#for plots:
plots<-c(1:29)

#the bit of code I used to figure out which species were missing in the 
#desig file.  Then I went in and manually added them and their information.
#a1<-levels(summerbymod$species) #or ''$genus
#a2<-levels(match$species) #or ''$genus
#different.names<- (!a1 %in% a2)
#not.in.a2<-a1[different.names]

#remove unused levels in desig$genus and desig$species
match<-droplevels(desig[desig$genus %in% levels(summerbymod$genus),])
match<-droplevels(match[match$species %in% levels(summerbymod$species),])

for(i in 1:dim(summerbymod)[1]){
  for(j in 1:dim(match)[1]){
    if(summerbymod$genus[i]==match$genus[j] && 
       summerbymod$species[i]==match$species[j])
    {summerbymod$form[i]<-as.character(match$form[j])}
  }
}

for(i in 1:dim(summerbymod)[1]){
  for(j in 1:dim(match)[1]){
    if(summerbymod$genus[i]==match$genus[j] && 
        summerbymod$species[i]==match$species[j])
    {summerbymod$month[i]<-as.character(match$designation[j])}
  }
}

#REMOVE SPRING FORBS?
if(designations == 1){
  summerspecies<-summerbymod[!summerbymod$month == "spring",]}else{
    summerspecies<-summerbymod
  }

#REMOVE NON-FORBS?
#summerspecies<-summerspecies[summerspecies$form == "forb",]

#RUN FROM HERE FOR SUBSEQUENT VERSIONS

otherforms<-0 # to remove forbs
forbsonly<-1 #1 to remove non-forbs
normalize<-0 #1 to return normalized (log- or sqrt-transformed data)
z<-0
corrstructure<-0 #1 to return visualizations

if(otherforms == 1){ 
  summerspecies<-summerspecies[summerspecies$form != "forb",]
}
#richness

#alpha
summmods.rich<-split(summerspecies[,c(7,9,11,13)],summerspecies[,1],drop=FALSE)
summalpha.rich<-matrix(data=NA,ncol=dim(summmods.rich[[1]])[2],nrow=length(summmods.rich))
for(i in 1:length(summmods.rich)){for(j in 1:dim(summmods.rich[[i]])[2]){
  summalpha.rich[i,j]<-sum(summmods.rich[[i]][j])
}
}
summalpha.rich[summalpha.rich == 0]<-NA
colnames(summalpha.rich)<-c("subplot1.rich","subplot2.rich","subplot3.rich",
                            "subplot4.rich")
summalpha.rich.keep<-summalpha.rich #keep a version in this format

#mean alpha
summmean.rich<-apply(summalpha.rich, 1, mean,na.rm=TRUE)

#gamma
summgamma.rich<-NULL
for(i in 1:length(summmods.rich)){summgamma.rich[i]<-dim(summmods.rich[[i]])[1]}
#normal

#beta
summbeta.rich<-summmean.rich/summgamma.rich
#normal

#variance
library(raster)
summvar.rich<-apply(summalpha.rich,1,cv,na.rm=TRUE)

#richness dataframe
#subplot
summalpha.rich<-as.matrix(unlist(as.list(t(summalpha.rich))))
summalpha.rich<-summalpha.rich[!is.na(summalpha.rich)]
#plot
summplot.rich<-data.frame(summgamma.rich,summmean.rich,
                          summbeta.rich,summvar.rich)

#NORMALIZE?
if(normalize==1){
  summgamma.rich.n<-sqrt(summgamma.rich)
  summmean.rich.n<-sqrt(summmean.rich)
  summvar.rich.n<-sqrt(summvar.rich)
  summplot.rich.n<-data.frame(summgamma.rich.n,
                              summmean.rich.n,summbeta.rich,
                              summvar.rich.n)
  summalpha.rich.n<-sqrt(summalpha.rich)
}

#evenness

midpoints<-c(.025,.05,1.5,3.5,7.5,17.5,37.5,62.5,85,97.5) #to convert cover classes
summmods.cover<-split(summerspecies[,c(6,8,10,12)],summerspecies[,1],drop=FALSE)

for(i in 1:length(summmods.cover)){
  for(j in 1:dim(summmods.cover[[i]])[2]){
    for(k in 1:dim(summmods.cover[[i]])[1]){if(is.na(summmods.cover[[i]][k,j])){summmods.cover[[i]][k,j]<-0}
      else(summmods.cover[[i]][k,j]<-midpoints[summmods.cover[[i]][k,j]])}
  }
}

summmods.cover.keep<-summmods.cover #keep a version of cover data for later

summmods.prop<-summmods.cover
for(i in 1:length(summmods.cover)){
  for(j in 1:dim(summmods.cover[[i]])[2]){
    for(k in 1:dim(summmods.cover[[i]])[1]){
      total<-sum(summmods.cover[[i]][,j])
      summmods.prop[[i]][k,j]<-(summmods.cover[[i]][k,j]/total)^2}
  }
}

summalpha.evens<-matrix(data=NA,ncol=dim(summmods.cover[[1]])[2],
                        nrow=length(summmods.cover),byrow=TRUE)
for(i in 1:length(summmods.cover)){for(j in 1:dim(summmods.cover[[i]])[2]){
  summalpha.evens[i,j]<-(1/sum(summmods.prop[[i]][,j]))
}
}

summalpha.evens<-(summalpha.evens/summalpha.rich.keep)
colnames(summalpha.evens)<-c("subplot1.evens","subplot2.evens","subplot3.evens",
                             "subplot4.evens")
summalpha.evens[summalpha.evens == 'NaN']<-NA

#mean alpha
summmean.evens<-apply(summalpha.evens, 1, mean,na.rm=TRUE)

#gamma
summplots.cover<-list()
for(i in 1:length(summmods.cover.keep)){
  summplots.cover[[i]]<-apply(summmods.cover.keep[[i]],1,sum)
}

summgamma.evens<-NULL
for(i in 1:length(summplots.cover)){
  total<-sum(summplots.cover[[i]])
  proportions<-NULL
  for(j in 1:length(summplots.cover[[i]])){
    proportions[j]<-(summplots.cover[[i]][j]/total)^2}
  summgamma.evens[i]<-1/sum(proportions)
}
summgamma.evens<-(summgamma.evens/summgamma.rich)

#beta
summbeta.evens<-(summmean.evens/summgamma.evens) #remember this won't be 0-1

#variance
summvar.evens<-apply(summalpha.evens,1,cv,na.rm=TRUE)

#evenness dataframes
#subplot
summalpha.evens<-as.matrix(unlist(as.list(t(summalpha.evens))))
summalpha.evens<-summalpha.evens[!is.na(summalpha.evens)]
#plot
summplot.evens<-data.frame(summgamma.evens,summmean.evens,
                           summbeta.evens,summvar.evens)

#NORMALIZE?
if(normalize==1){
  summmean.evens.n<-log(summmean.evens) #issues
  summgamma.evens.n<-log(summgamma.evens) 
  summbeta.evens.n<-log(summbeta.evens)
  summvar.evens.n<-sqrt(summvar.evens)
  summplot.evens.n<-data.frame(summgamma.evens.n,
                               summmean.evens.n,summbeta.evens.n,
                               summvar.evens.n)
  summalpha.evens.n<-log(summalpha.evens)
}

#all patterns
summplot.patterns<-data.frame(summplot.rich,summplot.evens)
summalpha.patterns<-data.frame(summalpha.rich,summalpha.evens)
if(z == 1){
  summplot.patterns.z<-scale(summplot.patterns, center=TRUE,scale=TRUE)
  summalpha.patterns.z<-scale(summalpha.patterns,center=TRUE,scale=TRUE)
}
#NORMALIZE?
if(normalize==1){
  summplot.patterns.n<-data.frame(summplot.rich.n,
                                  summplot.evens.n)
  summalpha.patterns.n<-data.frame(summalpha.rich.n,summalpha.evens.n)
  if(z == 1){
    summplot.patterns.n.z<-scale(summplot.patterns.n, center=TRUE,scale=TRUE)
    summalpha.patterns.n.z<-scale(summalpha.patterns.n,center=TRUE,scale=TRUE)
  }
}

#save versions:
if(z==0 & normalize == 0){
  summplot.patterns.otherforms<-summplot.patterns
  summalpha.patterns.otherforms<-summalpha.patterns
}else
  if(z==0 & normalize == 1){
    summplot.patterns.n.otherforms<-summplot.patterns.n
    summalpha.patterns.n.otherforms<-summalpha.patterns.n
  }else
    if(z==1 & normalize == 0){
      summplot.patterns.z.otherforms<-summplot.patterns.z
      summalpha.patterns.z.otherforms<-summalpha.patterns.z
    }else
      if(z==1 & normalize==1){
        summplot.patterns.n.z.otherforms<-summplot.patterns.n.z
        summalpha.patterns.n.z.otherforms<-summalpha.patterns.n.z
      }
#VISUALIZE CORRELATION STRUCTURE?
if(corrstructure==1){
  library(PerformanceAnalytics)
  chart.Correlation(summplot.rich.n)
  chart.Correlation(summplot.evens.n)
  chart.Correlation(summplot.patterns.n)
  chart.Correlation(summalpha.patterns.n)
}


if(forbsonly == 0){

#richness

#alpha
summmods.rich<-split(summerspecies[,c(7,9,11,13)],summerspecies[,1],drop=FALSE)
summalpha.rich<-matrix(data=NA,ncol=dim(summmods.rich[[1]])[2],nrow=length(summmods.rich))
for(i in 1:length(summmods.rich)){for(j in 1:dim(summmods.rich[[i]])[2]){
  summalpha.rich[i,j]<-sum(summmods.rich[[i]][j])
}
}
summalpha.rich[summalpha.rich == 0]<-NA
colnames(summalpha.rich)<-c("subplot1.rich","subplot2.rich","subplot3.rich",
                          "subplot4.rich")
summalpha.rich.keep<-summalpha.rich #keep a version in this format

#mean alpha
summmean.rich<-apply(summalpha.rich, 1, mean,na.rm=TRUE)

#gamma
summgamma.rich<-NULL
for(i in 1:length(summmods.rich)){summgamma.rich[i]<-dim(summmods.rich[[i]])[1]}
#normal

#beta
summbeta.rich<-summmean.rich/summgamma.rich
#normal

#variance
library(raster)
summvar.rich<-apply(summalpha.rich,1,cv,na.rm=TRUE)

#richness dataframe
#subplot
summalpha.rich<-as.matrix(unlist(as.list(t(summalpha.rich))))
summalpha.rich<-summalpha.rich[!is.na(summalpha.rich)]
#plot
summplot.rich<-data.frame(summgamma.rich,summmean.rich,
                              summbeta.rich,summvar.rich)

#NORMALIZE?
if(normalize==1){
  summgamma.rich.n<-sqrt(summgamma.rich)
  summmean.rich.n<-sqrt(summmean.rich)
  summvar.rich.n<-sqrt(summvar.rich)
  summplot.rich.n<-data.frame(summgamma.rich.n,
                                      summmean.rich.n,summbeta.rich,
                                      summvar.rich.n)
  summalpha.rich.n<-sqrt(summalpha.rich)
  }

#evenness

midpoints<-c(.025,.05,1.5,3.5,7.5,17.5,37.5,62.5,85,97.5) #to convert cover classes
summmods.cover<-split(summerspecies[,c(6,8,10,12)],summerspecies[,1],drop=FALSE)

for(i in 1:length(summmods.cover)){
  for(j in 1:dim(summmods.cover[[i]])[2]){
    for(k in 1:dim(summmods.cover[[i]])[1]){if(is.na(summmods.cover[[i]][k,j])){summmods.cover[[i]][k,j]<-0}
      else(summmods.cover[[i]][k,j]<-midpoints[summmods.cover[[i]][k,j]])}
  }
}

summmods.cover.keep<-summmods.cover #keep a version of cover data for later

summmods.prop<-summmods.cover
for(i in 1:length(summmods.cover)){
  for(j in 1:dim(summmods.cover[[i]])[2]){
    for(k in 1:dim(summmods.cover[[i]])[1]){
      total<-sum(summmods.cover[[i]][,j])
      summmods.prop[[i]][k,j]<-(summmods.cover[[i]][k,j]/total)^2}
  }
}

summalpha.evens<-matrix(data=NA,ncol=dim(summmods.cover[[1]])[2],
                      nrow=length(summmods.cover),byrow=TRUE)
for(i in 1:length(summmods.cover)){for(j in 1:dim(summmods.cover[[i]])[2]){
  summalpha.evens[i,j]<-(1/sum(summmods.prop[[i]][,j]))
  }
}

summalpha.evens<-(summalpha.evens/summalpha.rich.keep)
colnames(summalpha.evens)<-c("subplot1.evens","subplot2.evens","subplot3.evens",
                           "subplot4.evens")
summalpha.evens[summalpha.evens == 'NaN']<-NA

#mean alpha
summmean.evens<-apply(summalpha.evens, 1, mean,na.rm=TRUE)

#gamma
summplots.cover<-list()
for(i in 1:length(summmods.cover.keep)){
  summplots.cover[[i]]<-apply(summmods.cover.keep[[i]],1,sum)
}

summgamma.evens<-NULL
for(i in 1:length(summplots.cover)){
  total<-sum(summplots.cover[[i]])
  proportions<-NULL
  for(j in 1:length(summplots.cover[[i]])){
    proportions[j]<-(summplots.cover[[i]][j]/total)^2}
  summgamma.evens[i]<-1/sum(proportions)
}
summgamma.evens<-(summgamma.evens/summgamma.rich)

#beta
summbeta.evens<-(summmean.evens/summgamma.evens) #remember this won't be 0-1

#variance
summvar.evens<-apply(summalpha.evens,1,cv,na.rm=TRUE)

#evenness dataframes
#subplot
summalpha.evens<-as.matrix(unlist(as.list(t(summalpha.evens))))
summalpha.evens<-summalpha.evens[!is.na(summalpha.evens)]
#plot
summplot.evens<-data.frame(summgamma.evens,summmean.evens,
                               summbeta.evens,summvar.evens)

#NORMALIZE?
if(normalize==1){
  summmean.evens.n<-log(summmean.evens) #issues
  summgamma.evens.n<-log(summgamma.evens) 
  summbeta.evens.n<-log(summbeta.evens)
  summvar.evens.n<-sqrt(summvar.evens)
  summplot.evens.n<-data.frame(summgamma.evens.n,
                                      summmean.evens.n,summbeta.evens.n,
                                      summvar.evens.n)
  summalpha.evens.n<-log(summalpha.evens)
}

#all patterns
summplot.patterns<-data.frame(summplot.rich,summplot.evens)
summalpha.patterns<-data.frame(summalpha.rich,summalpha.evens)
if(z == 1){
  summplot.patterns.z<-scale(summplot.patterns, center=TRUE,scale=TRUE)
  summalpha.patterns.z<-scale(summalpha.patterns,center=TRUE,scale=TRUE)
}
#NORMALIZE?
if(normalize==1){
  summplot.patterns.n<-data.frame(summplot.rich.n,
                                      summplot.evens.n)
  summalpha.patterns.n<-data.frame(summalpha.rich.n,summalpha.evens.n)
  if(z == 1){
    summplot.patterns.n.z<-scale(summplot.patterns.n, center=TRUE,scale=TRUE)
    summalpha.patterns.n.z<-scale(summalpha.patterns.n,center=TRUE,scale=TRUE)
  }
}

#save versions:
if(z==0 & normalize == 0){
  summplot.patterns.allforms<-summplot.patterns
  summalpha.patterns.allforms<-summalpha.patterns
}else
  if(z==0 & normalize == 1){
    summplot.patterns.n.allforms<-summplot.patterns.n
    summalpha.patterns.n.allforms<-summalpha.patterns.n
  }else
    if(z==1 & normalize == 0){
      summplot.patterns.z.allforms<-summplot.patterns.z
      summalpha.patterns.z.allforms<-summalpha.patterns.z
    }else
      if(z==1 & normalize==1){
        summplot.patterns.n.z.allforms<-summplot.patterns.n.z
        summalpha.patterns.n.z.allforms<-summalpha.patterns.n.z
      }
#VISUALIZE CORRELATION STRUCTURE?
if(corrstructure==1){
  library(PerformanceAnalytics)
  chart.Correlation(summplot.rich.n)
  chart.Correlation(summplot.evens.n)
  chart.Correlation(summplot.patterns.n)
  chart.Correlation(summalpha.patterns.n)
  }
}

#REMOVE NON-FORBS?
if(forbsonly == 1){
  summerspecies<-summerspecies[summerspecies$form == "forb",]

#richness

#alpha
summmods.rich<-split(summerspecies[,c(7,9,11,13)],summerspecies[,1],drop=FALSE)
summalpha.rich<-matrix(data=NA,ncol=dim(summmods.rich[[1]])[2],nrow=length(summmods.rich))
for(i in 1:length(summmods.rich)){for(j in 1:dim(summmods.rich[[i]])[2]){
  summalpha.rich[i,j]<-sum(summmods.rich[[i]][j])
}
}
summalpha.rich[6,3:4]<-NA
colnames(summalpha.rich)<-c("subplot1.rich","subplot2.rich","subplot3.rich",
                            "subplot4.rich")
summalpha.rich.keep<-summalpha.rich #keep a version in this format

#mean alpha
summmean.rich<-apply(summalpha.rich, 1, mean,na.rm=TRUE)

#gamma
summgamma.rich<-NULL
for(i in 1:length(summmods.rich)){summgamma.rich[i]<-dim(summmods.rich[[i]])[1]}
#normal

#beta
summbeta.rich<-summmean.rich/summgamma.rich
#normal

#variance
library(raster)
summvar.rich<-apply(summalpha.rich,1,cv,na.rm=TRUE)

#richness dataframe
#subplot
summalpha.rich<-as.matrix(unlist(as.list(t(summalpha.rich))))
summalpha.rich<-summalpha.rich[!is.na(summalpha.rich)]

#plot
summplot.rich<-data.frame(summgamma.rich,summmean.rich,
                          summbeta.rich,summvar.rich)

#NORMALIZE?
if(normalize==1){
  summgamma.rich[summgamma.rich == 0]<-.000001
  summgamma.rich.n<-sqrt(summgamma.rich)
  summmean.rich[summmean.rich == 0]<-.000001
  summmean.rich.n<-sqrt(summmean.rich)
  summvar.rich[summvar.rich == 0]<-.000001
  summvar.rich.n<-sqrt(summvar.rich)
  summbeta.rich[summbeta.rich == 0]<-.000001
  summbeta.rich.n<-sqrt(summbeta.rich)
  summplot.rich.n<-data.frame(summgamma.rich.n,
                                 summmean.rich.n,summbeta.rich.n,
                                 summvar.rich.n)
  summalpha.rich[summalpha.rich == 0]<-.000001
  summalpha.rich.n<-sqrt(summalpha.rich)
  }

#evenness

midpoints<-c(.025,.05,1.5,3.5,7.5,17.5,37.5,62.5,85,97.5) #to convert cover classes
summmods.cover<-split(summerspecies[,c(6,8,10,12)],summerspecies[,1],drop=FALSE)

for(i in 1:length(summmods.cover)){
  for(j in 1:dim(summmods.cover[[i]])[2]){
    for(k in 1:dim(summmods.cover[[i]])[1]){if(is.na(summmods.cover[[i]][k,j])){summmods.cover[[i]][k,j]<-0}
      else(summmods.cover[[i]][k,j]<-midpoints[summmods.cover[[i]][k,j]])}
  }
}

summmods.cover.keep<-summmods.cover #keep a version of cover data for later

summmods.prop<-summmods.cover
for(i in 1:length(summmods.cover)){
  for(j in 1:dim(summmods.cover[[i]])[2]){
    for(k in 1:dim(summmods.cover[[i]])[1]){
      total<-sum(summmods.cover[[i]][,j])
      summmods.prop[[i]][k,j]<-(summmods.cover[[i]][k,j]/total)^2}
  }
}

summalpha.evens<-matrix(data=NA,ncol=dim(summmods.cover[[1]])[2],
                        nrow=length(summmods.cover),byrow=TRUE)
for(i in 1:length(summmods.cover)){for(j in 1:dim(summmods.cover[[i]])[2]){
  summalpha.evens[i,j]<-(1/sum(summmods.prop[[i]][,j]))
  }
}

summalpha.evens<-(summalpha.evens/summalpha.rich.keep)
colnames(summalpha.evens)<-c("subplot1.evens","subplot2.evens","subplot3.evens",
                             "subplot4.evens")

#mean alpha
summmean.evens<-apply(summalpha.evens, 1, mean,na.rm=TRUE)
#with richness I'd want to replace NaN values with 0 so that they're included in
#the mean.  But zero evenness doesn't make sense.  So they're excluded.

#gamma
summplots.cover<-list()
for(i in 1:length(summmods.cover.keep)){
  summplots.cover[[i]]<-apply(summmods.cover.keep[[i]],1,sum)
}

summgamma.evens<-NULL
for(i in 1:length(summplots.cover)){
  total<-sum(summplots.cover[[i]])
  proportions<-NULL
  for(j in 1:length(summplots.cover[[i]])){
    proportions[j]<-(summplots.cover[[i]][j]/total)^2}
  summgamma.evens[i]<-1/sum(proportions)
}
summgamma.evens<-(summgamma.evens/summgamma.rich)

#beta
summbeta.evens<-(summmean.evens/summgamma.evens)

#variance
summvar.evens<-apply(summalpha.evens,1,cv,na.rm=TRUE)

#evenness dataframes
#subplot
summalpha.evens<-as.matrix(unlist(as.list(t(summalpha.evens))))
summalpha.evens<-c(summalpha.evens[1:22],summalpha.evens[25:116])

#plot
summplot.evens<-data.frame(summgamma.evens,summmean.evens,
                           summbeta.evens,summvar.evens)

#NORMALIZE?
if(normalize==1){
  summmean.evens[summmean.evens == 0]<-.000001
  summmean.evens.n<-log(summmean.evens)
  summgamma.evens[summgamma.evens == 0]<-.000001
  summgamma.evens.n<-sqrt(summgamma.evens)
  summbeta.evens[summbeta.evens == 0]<-.000001 #can't salvage.
  summbeta.evens.n<-log(summbeta.evens)
  summvar.evens[summvar.evens == 0]<-.000001
  summvar.evens.n<-sqrt(summvar.evens)
  summplot.evens.n<-data.frame(summgamma.evens.n,
                                          summmean.evens.n,summbeta.evens,
                                          summvar.evens.n)
  summalpha.evens[summalpha.evens == 0]<-.000001
  summalpha.evens.n<-log(summalpha.evens)#can't really salvage this either.
}

#all patterns
summplot.patterns<-data.frame(summplot.rich,summplot.evens)
summalpha.patterns<-data.frame(summalpha.rich,summalpha.evens)
if(z == 1){
  summplot.patterns.z<-scale(summplot.patterns, center=TRUE,scale=TRUE)
  summalpha.patterns.z<-scale(summalpha.patterns,center=TRUE,scale=TRUE)
}
#NORMALIZE?
if(normalize==1){
  summplot.patterns.n<-data.frame(summplot.rich.n,
                                     summplot.evens.n)
  summalpha.patterns.n<-data.frame(summalpha.rich.n,summalpha.evens.n)
  if(z == 1){
    summplot.patterns.n.z<-scale(summplot.patterns.n, center=TRUE,scale=TRUE)
    summalpha.patterns.n.z<-scale(summalpha.patterns.n,center=TRUE,scale=TRUE)
  }
}
#save versions:
if(z==0 & normalize == 0){
  summplot.patterns.forbs<-summplot.patterns
  summalpha.patterns.forbs<-summalpha.patterns
}else
  if(z==0 & normalize == 1){
    summplot.patterns.n.forbs<-summplot.patterns.n
    summalpha.patterns.n.forbs<-summalpha.patterns.n
  }else
    if(z==1 & normalize == 0){
      summplot.patterns.z.forbs<-summplot.patterns.z
      summalpha.patterns.z.forbs<-summalpha.patterns.z
    }else
      if(z==1 & normalize==1){
        summplot.patterns.n.z.forbs<-summplot.patterns.n.z
        summalpha.patterns.n.z.forbs<-summalpha.patterns.n.z
      }

#VISUALIZE CORRELATION STRUCTURE?
if(corrstructure==1){
  library(PerformanceAnalytics)
  chart.Correlation(summplot.rich.n)
  chart.Correlation(summplot.evenness.n)
  chart.Correlation(summplot.patterns.n)
  chart.Correlation(summalpha.patterns.n)
  }
}