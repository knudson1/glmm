#check objfun using finite differences
library(glmm)
data(BoothHobert)
set.seed(1234)
out<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
family.glmm=bernoulli.glmm,m=50,doPQL=FALSE,debug=TRUE,distrib="normal")
mod.mcml<-out$mod.mcml
debug<-out$debug
nu.pql<-debug$nu.pql
nu.pql
nu.gen<-debug$nu.gen
beta.pql<-debug$beta.pql
beta.pql
family.glmm<-out$family.glmm
umat<-debug$umat
u.pql<-debug$u.star

par<-c(6,1.5)
del<-rep(.00000001,2)

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

# define a few things that will be used for finite differences
lth<-objfun(par=par, nbeta=1, nu.pql=nu.gen, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,distrib="normal",p1=p1,p2=p2,p3=p3, Sigmuh=Sigmuh,D.star=D.star)
lthdel<-objfun(par=par+del, nbeta=1, nu.pql=nu.pql, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,distrib="normal",p1=p1,p2=p2,p3=p3, Sigmuh=Sigmuh,D.star=D.star)

all.equal(as.vector(lth$gradient%*%del),lthdel$value-lth$value)
all.equal(as.vector(lth$hessian%*%del),lthdel$gradient-lth$gradient)




