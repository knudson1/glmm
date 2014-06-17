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

lth<-objfun(par=par, nbeta=1, nu.pql=nu.gen, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,distrib="normal")
lthdel<-objfun(par=par+del, nbeta=1, nu.pql=nu.pql, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,distrib="normal")

all.equal(as.vector(lth$gradient%*%del),lthdel$value-lth$value)
all.equal(as.vector(lth$hessian%*%del),lthdel$gradient-lth$gradient)




