#SPRING 2015 abiotic data
setwd("C:/Users/Sam/Desktop/PCAP")

#make some tough choices:
raw<-0 #do you want the raw values?
normalize<-1 #do you want to log or sqrt transform the data to fit a normal distribution?
z<-1 #do you want to transform the data to values between 0 and 1?

#soil
spsoil1<-read.csv("SPRING 2015 soil1.csv",header = TRUE,sep = ',')
spsoil2<-read.csv("SPRING 2015 soil2.csv",header = TRUE,sep = ',')

spsoil<-rbind(spsoil1,spsoil2)
spsoil<-spsoil[order((spsoil$plot.index),decreasing=FALSE),]

spalpha.soil<-spsoil[,2:16]
colnames(spalpha.soil)<-c("plot.index","mod","OM","P","Kppm","Mgppm","Cappm","pH","CEC","K","Mg",
                    "Ca","C","N","C.N")

#light
splight<-read.csv("SPRING 2015 light.csv",header = TRUE, sep = ',')
splightavg<-as.list(t(splight[,20:23]))
spalpha.light<-as.matrix(unlist(splightavg))

#litter, organic layer, and depth to restrictive layer
splor<-read.csv("SPRING 2015 litter_rd.csv",header = TRUE, sep = ',')
splalpha<-as.list(t(splor[,c(4,7,10,13)]))
spalpha.l<-as.matrix(unlist(splalpha))
sploalpha<-as.list(t(splor[,c(3,6,9,12)]))
spalpha.lo<-as.matrix(unlist(sploalpha))
spalpha.o<-(spalpha.lo-spalpha.l)
spralpha<-as.list(t(splor[,c(5,8,11,14)]))
spalpha.r<-as.matrix(unlist(spralpha))
spalpha.lor<-data.frame(spalpha.l,spalpha.o,spalpha.r)

#RAW?
if(raw == 1){
#subplot
spalpha.ab<-data.frame(spalpha.soil,I(na.omit(spalpha.light)),I(na.omit(spalpha.l)),I(na.omit(spalpha.o)),
                       I(na.omit(spalpha.r)))
colnames(spalpha.ab)<-c(names(spalpha.soil),"spalpha.light","spalpha.l","spalpha.o","spalpha.r")

#plot
#mean & variance
library(raster)
spabbymod<-split(spalpha.ab[,3:19],spalpha.ab$plot.index,drop=FALSE)
spmean.ab<-matrix(data=NA,nrow=length(spabbymod),ncol=dim(spabbymod[[1]])[2],byrow = TRUE)
spvar.ab<-matrix(data=NA,nrow=length(spabbymod),ncol=dim(spabbymod[[1]])[2],byrow = TRUE)
for(i in 1:length(spabbymod)){
  for(j in 1:dim(spabbymod[[i]])[2]){
    spmean.ab[i,j]<-mean(as.numeric(unlist(spabbymod[[i]][j])))
    spvar.ab[i,j]<-cv(as.numeric(unlist(spabbymod[[i]][j])))
    if(is.na(spvar.ab[i,j])){spvar.ab[i,j]<-0}}
}

colnames(spmean.ab)<-c('mOM','mP','mKppm','mMgppm','mCappm','mpH','mCEC',
                         'mK','mMg','mCa','mC','mN','mC.N','mlight','ml',
                         'mo','mr')
colnames(spvar.ab)<-c('vOM','vP','vKppm','vMgppm','vCappm','vpH','vCEC',
                        'vK','vMg','vCa','vC','vN','vC.N','vlight','vl',
                        'vo','vr')

spplot.ab<-data.frame(spmean.ab,spvar.ab)

if(z==1){
  plots.mods<-spalpha.ab[,1:2] #mods by plots to tack back on later
  temp<-spalpha.ab[,c(-1,-2)]
  temp.z<-list()
  for(i in 1:dim(temp)[2]){
    temp.z[[i]]<-scale(temp[,i],center=TRUE,scale=TRUE)
  }
  spalpha.ab<-data.frame(matrix(unlist(temp.z),nrow=dim(temp)[1],ncol=dim(temp)[2]))
  colnames(spalpha.ab)<-c(names(temp))
  
  spmean.ab.z<-list()
  spvar.ab.z<-list()
  for(i in 1:dim(spmean.ab)[2]){
    spmean.ab.z[[i]]<-scale(spmean.ab[,i],scale=TRUE,center=TRUE)
    spvar.ab.z[[i]]<-scale(spvar.ab[,i],scale=TRUE,center=TRUE)
  }
  spmean.ab.z<-data.frame(matrix(unlist(spmean.ab.z),nrow=dim(spmean.ab)[1],ncol=dim(spmean.ab)[2]))
  colnames(spmean.ab.z)<-c('mOM','mP','mKppm','mMgppm','mCappm','mpH','mCEC',
                               'mK','mMg','mCa','mC','mN','mC.N','mlight','ml',
                               'mo','mr')
  spvar.ab.z<-data.frame(matrix(unlist(spvar.ab.z),nrow=dim(spvar.ab)[1],ncol=dim(spvar.ab)[2]))
  colnames(spvar.ab.z)<-c('vOM','vP','vKppm','vMgppm','vCappm','vpH','vCEC',
                              'vK','vMg','vCa','vC','vN','vC.N','vlight','vl',
                              'vo','vr')
  
  spplot.ab.z<-data.frame(spmean.ab.z,spvar.ab.z)
  spalpha.ab.z<-data.frame(plots.mods,spalpha.ab)
}
}

#NORMALIZE?
if(normalize == 1){
  #subplot
  spalpha.ab<-data.frame(spalpha.soil,I(na.omit(spalpha.light)),I(na.omit(spalpha.l)),I(na.omit(spalpha.o)),
                         I(na.omit(spalpha.r)))
  colnames(spalpha.ab)<-c(names(spalpha.soil),"spalpha.light","spalpha.l","spalpha.o","spalpha.r")
  
  #need to replace 0s to log transform:
  spalpha.ab[spalpha.ab == 0]<-.000001
  temp<-spalpha.ab[,3:19]
  for(i in 1:dim(temp)[2]){
    if(i==2|i==5|i==6|i==8|i==13|i==16){temp[,i]<-log(temp[,i])}
      else if(i==1|i==7|i==9){temp[,i]<-temp[,i]}else{
        temp[,i]<-sqrt(temp[,i])}
    }
  spalpha.ab.t<-data.frame(spalpha.ab[,1:2],temp)

  #plot
  #mean & variance
  library(raster)
  spabbymod<-split(spalpha.ab[,3:19],spalpha.ab$plot.index,drop=FALSE)
  spmean.ab<-matrix(data=NA,nrow=length(spabbymod),ncol=dim(spabbymod[[1]])[2],byrow = TRUE)
  spvar.ab<-matrix(data=NA,nrow=length(spabbymod),ncol=dim(spabbymod[[1]])[2],byrow = TRUE)
  for(i in 1:length(spabbymod)){
    for(j in 1:dim(spabbymod[[i]])[2]){
      spmean.ab[i,j]<-mean(as.numeric(unlist(spabbymod[[i]][j])))
      spvar.ab[i,j]<-cv(as.numeric(unlist(spabbymod[[i]][j])))
      if(is.na(spvar.ab[i,j])){spvar.ab[i,j]<-0}}
  }
  
  colnames(spmean.ab)<-c('mOM','mP','mKppm','mMgppm','mCappm','mpH','mCEC',
                         'mK','mMg','mCa','mC','mN','mC.N','mlight','ml',
                         'mo','mr')
  colnames(spvar.ab)<-c('vOM','vP','vKppm','vMgppm','vCappm','vpH','vCEC',
                        'vK','vMg','vCa','vC','vN','vC.N','vlight','vl',
                        'vo','vr')
  
  #transform mean
  #need to replace 0s to log transform:
  spmean.ab[spmean.ab == 0]<-.000001
  temp<-spmean.ab
  for(i in 1:dim(temp)[2]){
    if(i==2|i==5|i==6|i==8|i==13|i==16){temp[,i]<-log(temp[,i])}
    else if(i==1|i==7|i==9){temp[,i]<-temp[,i]}else{
      temp[,i]<-sqrt(temp[,i])}}
  spmean.ab.t<-temp
  
#transform variance
  spvar.ab[spvar.ab == 0]<-.000001
  temp<-spvar.ab
  for(i in 1:dim(temp)[2]){
    if(i==2|i==3|i==5|i==9|i==10|i==16|i==17){temp[,i]<-sqrt(temp[,i])}
    else if(i==1){temp[,i]<-temp[,i]}else{
      temp[,i]<-log(temp[,i])}}
  spvar.ab.t<-temp
  
  spplot.ab.t<-data.frame(spmean.ab.t,spvar.ab.t)

  if(z==1){
    plots.mods<-spalpha.ab.t[,1:2]
    temp<-spalpha.ab.t[,c(-1,-2)]
    temp.z<-list()
    for(i in 1:dim(temp)[2]){
      temp.z[[i]]<-scale(temp[,i],center=TRUE,scale=TRUE)
    }
    spalpha.ab.t.z<-data.frame(matrix(unlist(temp.z),nrow=dim(temp)[1],ncol=dim(temp)[2]))
    colnames(spalpha.ab.t.z)<-c(names(temp))

    spmean.ab.t.z<-list()
    spvar.ab.t.z<-list()
    for(i in 1:dim(spmean.ab.t)[2]){
      spmean.ab.t.z[[i]]<-scale(spmean.ab.t[,i],scale=TRUE,center=TRUE)
      spvar.ab.t.z[[i]]<-scale(spvar.ab.t[,i],scale=TRUE,center=TRUE)
    }
    spmean.ab.t.z<-data.frame(matrix(unlist(spmean.ab.t.z),nrow=dim(spmean.ab.t)[1],ncol=dim(spmean.ab.t)[2]))
    colnames(spmean.ab.t.z)<-c('mOM','mP','mKppm','mMgppm','mCappm','mpH','mCEC',
                           'mK','mMg','mCa','mC','mN','mC.N','mlight','ml',
                           'mo','mr')
    spvar.ab.t.z<-data.frame(matrix(unlist(spvar.ab.t.z),nrow=dim(spvar.ab.t)[1],ncol=dim(spvar.ab.t)[2]))
    colnames(spvar.ab.t.z)<-c('vOM','vP','vKppm','vMgppm','vCappm','vpH','vCEC',
                          'vK','vMg','vCa','vC','vN','vC.N','vlight','vl',
                          'vo','vr')
    
    spplot.ab.t.z<-data.frame(spmean.ab.t.z,spvar.ab.t.z)
    spalpha.ab.t.z<-data.frame(plots.mods,spalpha.ab.t.z)
  }
}

#delta.light
setwd("C:/Users/Sam/Desktop/PCAP")
splight<-read.csv("SPRING 2015 light.csv",header = TRUE, sep = ',')
splight.delta<-splight[,3:18]

#to get one variance for the entire plot, 