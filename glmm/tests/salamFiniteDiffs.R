library(glmm)
set.seed(1234)
data(salamander)
sal<-glmm(Mate~Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"), data=salamander,family.glmm=bernoulli.glmm,m=100,debug=TRUE,doPQL=FALSE, cores=2)

vars <- new.env(parent = emptyenv())
objfun<-glmm:::objfun
beta<-rep(0,4)
nu<-rep(2,2)
par<-c(beta,nu)
vars$nbeta<-4

debug<-sal$debug
vars$nu.pql<-debug$nu.pql
beta.pql<-debug$beta.pql
vars$umat<-debug$umat
vars$m1<-debug$m1
vars$p1<-vars$p2<-vars$p3<-1/3
vars$mod.mcml<-sal$mod.mcml
vars$u.star<-debug$u.star
vars$family.glmm=bernoulli.glmm
vars$zeta<-5

getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

#need to get D*
eek<-getEk(vars$mod.mcml$z)
Aks<-Map("*",eek,vars$nu.pql)
vars$D.star<-addVecs(Aks) 
vars$D.star<-diag(vars$D.star)
D.star.inv<-solve(vars$D.star)

cache<-new.env(parent = emptyenv())

#need to also recreate the variance matrix of last imp sampling distribution
Z=do.call(cbind,vars$mod.mcml$z)
eta.star<-as.vector(vars$mod.mcml$x%*%beta.pql+Z%*%vars$u.star)
cdouble<-bernoulli.glmm()$cpp(eta.star) #still a vector
cdouble<-diag(cdouble)
vars$Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
vars$Sigmuh<-solve(vars$Sigmuh.inv)

#evaluate objective function at two close places
del<-rep(10^-6,6)
vars$ntrials<-1

vars$no_cores <- sal$cores

ltheta<-objfun(par, cache, vars=vars)

lthetadel<-objfun(par+del, cache, vars=vars)

#do finite diffs to check value --> gradient
c(as.vector(ltheta$gradient%*%del),lthetadel$value-ltheta$value)
as.vector(ltheta$gradient%*%del)-(lthetadel$value-ltheta$value)
all.equal(as.vector(ltheta$gradient%*%del),lthetadel$value-ltheta$value,tolerance=10^-5)

#do finite diffs to check gradient --> hessian 
all.equal(lthetadel$gradient-ltheta$gradient,as.vector(ltheta$hessian%*%del),tolerance=10^-6)
cbind(lthetadel$gradient-ltheta$gradient,as.vector(ltheta$hessian%*%del))
#how big are the diffs? very small
lthetadel$gradient-ltheta$gradient-as.vector(ltheta$hessian%*%del)


