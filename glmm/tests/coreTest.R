library(glmm)
data(BoothHobert)
set.seed(1234)
out<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
          family.glmm=bernoulli.glmm,m=50,doPQL=FALSE,debug=TRUE, cores=2)
vars <- new.env(parent = emptyenv())
vars$mod.mcml<-out$mod.mcml
debug<-out$debug
vars$nu.pql<-debug$nu.pql
vars$nu.pql
beta.pql<-debug$beta.pql
beta.pql
vars$family.glmm<-out$family.glmm
vars$umat<-debug$umat
u.pql<-debug$u.star
vars$m1<-debug$m1
vars$ntrials<-1

par<-c(6,1.5)
del<-rep(10^-8,2)

objfun<-glmm:::objfun
getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

#need to get D*
eek<-getEk(vars$mod.mcml$z)
Aks<-Map("*",eek,vars$nu.pql)
vars$D.star<-addVecs(Aks) 
vars$D.star<-diag(vars$D.star)
D.star.inv<-solve(vars$D.star)


#need to also recreate the variance matrix of last imp sampling distribution
Z=do.call(cbind,vars$mod.mcml$z)
eta.star<-as.vector(vars$mod.mcml$x%*%beta.pql+Z%*%u.pql)
cdouble<-bernoulli.glmm()$cpp(eta.star) #still a vector
cdouble<-diag(cdouble)
vars$Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
vars$Sigmuh<-solve(vars$Sigmuh.inv)

vars$p1=vars$p2=vars$p3=1/3
vars$zeta=5

vars$no_cores <- 1

vars$nbeta <- 1
vars$u.star <- u.pql

core1<-objfun(par=par, vars=vars)

vars$no_cores <- out$cores

core2<-objfun(par=par, vars=vars)

all.equal(core1, core2)
