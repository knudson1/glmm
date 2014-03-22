# test on Booth and Hobert data

source("MCLA.R")
dat<-read.csv("BoothHobertData.csv")
 dat$z1<-as.factor(dat$z1)

#mostly just want to get the PQL estimate 
mod.mcml<-mcml(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=dat,family.mcml=bernoulli.mcml,m=4) 