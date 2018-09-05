library(glmm)
data(BoothHobert)
set.seed(1234)
out<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
          family.glmm=bernoulli.glmm,m=50,doPQL=FALSE,debug=TRUE, cores=1)
mod.mcml<-out$mod.mcml
debug<-out$debug
nu.pql<-debug$nu.pql
nu.pql
beta.pql<-debug$beta.pql
beta.pql
family.glmm<-out$family.glmm
umat<-debug$umat
u.pql<-debug$u.star
m1<-debug$m1
ntrials<-1

par<-c(6,1.5)
del<-rep(10^-8,2)

objfun<-glmm:::objfun
getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

#need to get D*
eek<-getEk(mod.mcml$z)
Aks<-Map("*",eek,nu.pql)
D.star<-addVecs(Aks) 
D.star<-diag(D.star)
D.star.inv<-solve(D.star)


#need to also recreate the variance matrix of last imp sampling distribution
Z=do.call(cbind,mod.mcml$z)
eta.star<-as.vector(mod.mcml$x%*%beta.pql+Z%*%u.pql)
cdouble<-bernoulli.glmm()$cpp(eta.star) #still a vector
cdouble<-diag(cdouble)
Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
Sigmuh<-solve(Sigmuh.inv)

p1=p2=p3=1/3
zeta=5

no_cores <- out$cores

core1<-objfun(par=par, nbeta=1, nu.pql=nu.pql, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,p1=p1,p2=p2,p3=p3,m1=m1, Sigmuh=Sigmuh, D.star=D.star, Sigmuh.inv= Sigmuh.inv, zeta=zeta, ntrials=ntrials, no_cores=1)
core2<-objfun(par=par, nbeta=1, nu.pql=nu.pql, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,p1=p1,p2=p2,p3=p3,m1=m1, Sigmuh=Sigmuh, D.star=D.star, Sigmuh.inv= Sigmuh.inv, zeta=zeta, ntrials=ntrials, no_cores=no_cores)

all.equal(core1, core2)