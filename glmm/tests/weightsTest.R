library(glmm)
data(BoothHobert)

set.seed(1234)
zeroweight<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
          family.glmm=bernoulli.glmm,m=50,doPQL=FALSE,debug=TRUE,weights=c(rep(0,3),rep(1,147)))

BoothHobert2 <- BoothHobert[-c(1:3),]
set.seed(1234)
removed<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert2,
                 family.glmm=bernoulli.glmm,m=50,doPQL=FALSE,debug=TRUE)

all.equal(coef(zeroweight),coef(removed))
all.equal(zeroweight$likelihood.value,removed$likelihood.value)
all.equal(zeroweight$likelihood.gradient,removed$likelihood.gradient)
all.equal(zeroweight$likelihood.hessian,removed$likelihood.hessian)

set.seed(1234)
doubleweight<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
                 family.glmm=bernoulli.glmm,m=50,doPQL=FALSE,debug=TRUE,weights=c(rep(1,149),2))

#umatadded <- doubleweight$umat

#nums <- c(rep(1,150),2)
#for(i in 1:nrow(BoothHobert)){
  #BoothHobert3[i] <- BoothHobert[i,]*nums[i]
#}
BoothHobert3 <- rbind(BoothHobert, BoothHobert[150,])
set.seed(1234)
added<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert3,
              family.glmm=bernoulli.glmm,m=50,doPQL=FALSE,debug=TRUE)

all.equal(coef(doubleweight),coef(added))
all.equal(doubleweight$likelihood.value,added$likelihood.value)
all.equal(doubleweight$likelihood.gradient,added$likelihood.gradient)
all.equal(doubleweight$likelihood.hessian,added$likelihood.hessian)