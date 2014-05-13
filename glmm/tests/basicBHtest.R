library(glmm)
data(BoothHobert)
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,family.glmm=bernoulli.glmm,m=100,doPQL=TRUE)
summary(mod.mcml1)
