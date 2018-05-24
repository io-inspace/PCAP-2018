#SUMMER 2015 abiotic data
setwd("C:/Users/Sam/Desktop/PCAP")

#make some tough choices:
raw<-0 #do you want the raw values?
normalize<-1 #do you want to log or sqrt transform the data to fit a normal distribution?
z<-1 #do you want to transform the data to values between 0 and 1?

#soil
summsoil1<-read.csv("SUMMER 2015 soil1.csv",header = TRUE,sep = ',')
summsoil2<-read.csv("SUMMER 2015 soil2.csv",header = TRUE,sep = ',')

summsoil<-rbind(summsoil1,summsoil2)
summsoil<-summsoil[order((summsoil$plot.index),decreasing=FALSE),]

summalpha.soil<-summsoil[,2:16]
colnames(summalpha.soil)<-c("plot.index","mod","OM","P","Kppm","Mgppm","Cappm","pH","CEC","K","Mg",
                    "Ca","C","N","C.N")

#there were two samples from plot 16:
double<-summalpha.soil[summalpha.soil$plot.index == 16,]
doubleordered<-double[order(double$mod,decreasing=FALSE),]
doublesplit<-split(doubleordered[,3:15],doubleordered$mod,drop=FALSE)
sixteen<-matrix(data=NA,ncol=dim(doublesplit[[1]])[2],nrow=length(doublesplit))
for(i in 1:length(doublesplit)){
  for(j in 1:dim(doublesplit[[i]])[2]){
    sixteen[i,j]<-mean(unlist(doublesplit[[i]][j]))
    }
}
sixteen<-data.frame(rep(16,4),c('a','b','c','d'),sixteen)
colnames(sixteen)<-c(names(summalpha.soil))
summalpha.soil<-rbind(summalpha.soil[!summalpha.soil$plot.index == 16,],sixteen)
summalpha.soil<-summalpha.soil[order(summalpha.soil$plot.index,decreasing = FALSE),]

#light
summlight<-read.csv("SUMMER 2015 light.csv",header = TRUE,sep = ',')
summlightavg<-as.list(t(summlight[,20:23]))
summalpha.light<-as.matrix(unlist(summlightavg))

#litter, organic layer, and depth to restrictive layer
summlor<-read.csv("SUMMER 2015 litter_rd.csv",header = TRUE, sep = ',')
summlalpha<-as.list(t(summlor[,c(4,7,10,13)]))
summalpha.l<-as.matrix(unlist(summlalpha))
summloalpha<-as.list(t(summlor[,c(3,6,9,12)]))
summalpha.lo<-as.matrix(unlist(summloalpha))
summalpha.o<-(summalpha.lo-summalpha.l)
summralpha<-as.list(t(summlor[,c(5,8,11,14)]))
summalpha.r<-as.matrix(unlist(summralpha))
summalpha.lor<-data.frame(summalpha.l,summalpha.o,summalpha.r)

#RAW?
if(raw == 1){
  #subplot
  summalpha.ab<-data.frame(summalpha.soil,I(na.omit(summalpha.light)),I(na.omit(summalpha.l)),I(na.omit(summalpha.o)),
                         I(na.omit(summalpha.r)))
  colnames(summalpha.ab)<-c(names(summalpha.soil),"summalpha.light","summalpha.l","summalpha.o","summalpha.r")
  
  #plot
  #mean & variance
  library(raster)
  summabbymod<-split(summalpha.ab[,3:19],summalpha.ab$plot.index,drop=FALSE)
  summmean.ab<-matrix(data=NA,nrow=length(summabbymod),ncol=dim(summabbymod[[1]])[2],byrow = TRUE)
  summvar.ab<-matrix(data=NA,nrow=length(summabbymod),ncol=dim(summabbymod[[1]])[2],byrow = TRUE)
  for(i in 1:length(summabbymod)){
    for(j in 1:dim(summabbymod[[i]])[2]){
      summmean.ab[i,j]<-mean(as.numeric(unlist(summabbymod[[i]][j])))
      summvar.ab[i,j]<-cv(as.numeric(unlist(summabbymod[[i]][j])))
      if(is.na(summvar.ab[i,j])){summvar.ab[i,j]<-0}}
  }
  
  colnames(summmean.ab)<-c('mOM','mP','mKppm','mMgppm','mCappm','mpH','mCEC',
                         'mK','mMg','mCa','mC','mN','mC.N','mlight','ml',
                         'mo','mr')
  colnames(summvar.ab)<-c('vOM','vP','vKppm','vMgppm','vCappm','vpH','vCEC',
                        'vK','vMg','vCa','vC','vN','vC.N','vlight','vl',
                        'vo','vr')
  
  summplot.ab<-data.frame(summmean.ab,summvar.ab)
  
  if(z==1){
    plots.mods<-summalpha.ab[,1:2] #mods by plots to tack back on later
    temp<-summalpha.ab[,c(-1,-2)]
    temp.z<-list()
    for(i in 1:dim(temp)[2]){
      temp.z[[i]]<-scale(temp[,i],center=TRUE,scale=TRUE)
    }
    summalpha.ab<-data.frame(matrix(unlist(temp.z),nrow=dim(temp)[1],ncol=dim(temp)[2]))
    colnames(summalpha.ab)<-c(names(temp))
    
    summmean.ab.z<-list()
    summvar.ab.z<-list()
    for(i in 1:dim(summmean.ab)[2]){
      summmean.ab.z[[i]]<-scale(summmean.ab[,i],scale=TRUE,center=TRUE)
      summvar.ab.z[[i]]<-scale(summvar.ab[,i],scale=TRUE,center=TRUE)
    }
    summmean.ab.z<-data.frame(matrix(unlist(summmean.ab.z),nrow=dim(summmean.ab)[1],ncol=dim(summmean.ab)[2]))
    colnames(summmean.ab.z)<-c('mOM','mP','mKppm','mMgppm','mCappm','mpH','mCEC',
                             'mK','mMg','mCa','mC','mN','mC.N','mlight','ml',
                             'mo','mr')
    summvar.ab.z<-data.frame(matrix(unlist(summvar.ab.z),nrow=dim(summvar.ab)[1],ncol=dim(summvar.ab)[2]))
    colnames(summvar.ab.z)<-c('vOM','vP','vKppm','vMgppm','vCappm','vpH','vCEC',
                            'vK','vMg','vCa','vC','vN','vC.N','vlight','vl',
                            'vo','vr')
    
    summplot.ab.z<-data.frame(summmean.ab.z,summvar.ab.z)
    summalpha.ab.z<-data.frame(plots.mods,summalpha.ab)
  }
}

#NORMALIZE?
if(normalize == 1){
  #subplot
  summalpha.ab<-data.frame(summalpha.soil,I(na.omit(summalpha.light)),I(na.omit(summalpha.l)),I(na.omit(summalpha.o)),
                         I(na.omit(summalpha.r)))
  colnames(summalpha.ab)<-c(names(summalpha.soil),"summalpha.light","summalpha.l","summalpha.o","summalpha.r")
  
  #need to replace 0s to log transform:
  summalpha.ab[summalpha.ab == 0]<-.000001
  temp<-summalpha.ab[,3:19]
  for(i in 1:dim(temp)[2]){
    if(i==3|i==9|i==10|i==14|i==15|i==17){temp[,i]<-sqrt(temp[,i])}
    else if(i==16){temp[,i]<-temp[,i]}
    else{
      temp[,i]<-log(temp[,i])}
  }
  summalpha.ab.t<-data.frame(summalpha.ab[,1:2],temp)
  #organic layer depth is so zero inflated it can't be transformed and 
  #shouldn't be included in most analyses, but I'm leaving it for 
  #comparison to spring.

  #plot
  #mean & variance
  library(raster)
  summabbymod<-split(summalpha.ab[,3:19],summalpha.ab$plot.index,drop=FALSE)
  summmean.ab<-matrix(data=NA,nrow=length(summabbymod),ncol=dim(summabbymod[[1]])[2],byrow = TRUE)
  summvar.ab<-matrix(data=NA,nrow=length(summabbymod),ncol=dim(summabbymod[[1]])[2],byrow = TRUE)
  for(i in 1:length(summabbymod)){
    for(j in 1:dim(summabbymod[[i]])[2]){
      summmean.ab[i,j]<-mean(as.numeric(unlist(summabbymod[[i]][j])))
      summvar.ab[i,j]<-cv(as.numeric(unlist(summabbymod[[i]][j])))
      if(is.na(summvar.ab[i,j])){summvar.ab[i,j]<-0}}
  }
  
  colnames(summmean.ab)<-c('mOM','mP','mKppm','mMgppm','mCappm','mpH','mCEC',
                         'mK','mMg','mCa','mC','mN','mC.N','mlight','ml',
                         'mo','mr')
  colnames(summvar.ab)<-c('vOM','vP','vKppm','vMgppm','vCappm','vpH','vCEC',
                        'vK','vMg','vCa','vC','vN','vC.N','vlight','vl',
                        'vo','vr')
  
  #transform mean
  #need to replace 0s to log transform:
  summmean.ab[summmean.ab == 0]<-.000001
  temp<-summmean.ab
  for(i in 1:dim(temp)[2]){
    if(i==3|i==14|i==15|i==17){temp[,i]<-sqrt(temp[,i])}
    else if(i==7|i==9|i==16){temp[,i]<-temp[,i]}else{
      temp[,i]<-log(temp[,i])}}
  summmean.ab.t<-temp
  
  #transform variance
  summvar.ab[summvar.ab == 0]<-.000001
  temp<-summvar.ab
  for(i in 1:dim(temp)[2]){
    if(i==2|i==3|i==8|i==13){temp[,i]<-log(temp[,i])}
    else if(i==16){temp[,i]<-temp[,i]}else{
      temp[,i]<-sqrt(temp[,i])}}
  summvar.ab.t<-temp
  #organic layer and rd can't be transformed and probably shouldn't 
  #be used

  summplot.ab.t<-data.frame(summmean.ab.t,summvar.ab.t)
  
  if(z==1){
    plots.mods<-summalpha.ab.t[,1:2]
    temp<-summalpha.ab.t[,c(-1,-2)]
    temp.z<-list()
    for(i in 1:dim(temp)[2]){
      temp.z[[i]]<-scale(temp[,i],center=TRUE,scale=TRUE)
    }
    summalpha.ab.t.z<-data.frame(matrix(unlist(temp.z),nrow=dim(temp)[1],ncol=dim(temp)[2]))
    colnames(summalpha.ab.t.z)<-c(names(temp))
    
    summmean.ab.t.z<-list()
    summvar.ab.t.z<-list()
    for(i in 1:dim(summmean.ab.t)[2]){
      summmean.ab.t.z[[i]]<-scale(summmean.ab.t[,i],scale=TRUE,center=TRUE)
      summvar.ab.t.z[[i]]<-scale(summvar.ab.t[,i],scale=TRUE,center=TRUE)
    }
    summmean.ab.t.z<-data.frame(matrix(unlist(summmean.ab.t.z),nrow=dim(summmean.ab.t)[1],ncol=dim(summmean.ab.t)[2]))
    colnames(summmean.ab.t.z)<-c('mOM','mP','mKppm','mMgppm','mCappm','mpH','mCEC',
                               'mK','mMg','mCa','mC','mN','mC.N','mlight','ml',
                               'mo','mr')
    summvar.ab.t.z<-data.frame(matrix(unlist(summvar.ab.t.z),nrow=dim(summvar.ab.t)[1],ncol=dim(summvar.ab.t)[2]))
    colnames(summvar.ab.t.z)<-c('vOM','vP','vKppm','vMgppm','vCappm','vpH','vCEC',
                              'vK','vMg','vCa','vC','vN','vC.N','vlight','vl',
                              'vo','vr')
    
    summplot.ab.t.z<-data.frame(summmean.ab.t.z,summvar.ab.t.z)
    summalpha.ab.t.z<-data.frame(plots.mods,summalpha.ab.t.z)
  }
}

