#2015 CATEGORICAL DATA

setwd("C:/Users/Sam/Desktop/PCAP")

categorical<-read.csv('2015 categorical data.csv', header = TRUE, sep = ',')

table(categorical[,3:4])