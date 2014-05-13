#check objfun using finite differences
library(glmm)
data(BoothHobert)
set.seed(1234)
out<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
family.glmm=bernoulli.glmm,m=1000,doPQL=FALSE)
mod.mcml<-out$mod.mcml
nu.pql<-out$nu.pql
beta.pql<-out$beta.pql
family.glmm<-out$family.glmm
umat<-out$umat
u.pql<-out$u.pql

par<-c(6,1.5)
del<-rep(.00000001,2)

objfun<-glmm:::objfun

lth<-objfun(par=par,nbeta=1,nu.pql=nu.pql,umat=umat,u.star=u.pql,mod.mcml=mod.mcml,family.glmm=family.glmm)
lthdel<-objfun(par=par+del,nbeta=1,nu.pql=nu.pql,umat=umat,u.star=u.pql,mod.mcml=mod.mcml,family.glmm=family.glmm)

all.equal(as.vector(lth$gradient%*%del),lthdel$value-lth$value)
all.equal(as.vector(lth$hessian%*%del),lthdel$gradient-lth$gradient)


