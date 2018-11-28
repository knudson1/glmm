library(glmm)
data(BoothHobert)
clust <- makeCluster(2)
set.seed(1234)
out<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
          family.glmm=bernoulli.glmm,m=50,doPQL=FALSE,debug=TRUE, cluster=clust)

vars <- new.env(parent = emptyenv())
debug<-out$debug

vars$m1 <- debug$m1
m2 <- debug$m2
m3 <- debug$m3

vars$zeta <- 5

vars$cl <- clust
registerDoParallel(vars$cl)                   #making cluster usable with foreach
vars$no_cores <- length(vars$cl)

vars$mod.mcml<-out$mod.mcml

vars$nu.pql <- debug$nu.pql

getEk<-glmm:::getEk
addVecs<-glmm:::addVecs
genRand<-glmm:::genRand

nrand<-lapply(vars$mod.mcml$z,ncol)
nrandom<-unlist(nrand)
q<-sum(nrandom)
totnrandom<-sum(nrandom)
s.pql<-rep(0,totnrandom)
if(q!=length(s.pql)) stop("Can't happen. Number of random effects returned by PQL must match number of random effects specified by model.")
eek<-getEk(vars$mod.mcml$z)
sigma.pql<-rep(1,length(vars$mod.mcml$z))
#if any of the variance components are too close to 0, make them bigger:
if(any(sigma.pql<10^-3)){
  theseguys<-which(sigma.pql<10^-3)
  sigma.pql[theseguys]<-10^-3
}
Aks<-Map("*",eek,sigma.pql)
A.star<-addVecs(Aks) #at this point still a vector
vars$D.star<-A.star*A.star #still a vector
vars$u.star<-A.star*s.pql 
Dstarinvdiag<-1/vars$D.star
Dstarnotsparse<-diag(vars$D.star)
D.star.inv<-Diagonal(length(vars$u.star),Dstarinvdiag)
vars$D.star<-Diagonal(length(vars$u.star),vars$D.star)

vars$family.glmm<-out$family.glmm

vars$ntrials<-1

beta.pql <- debug$beta.pql

simulate <- function(vars, Dstarnotsparse, m2, m3, beta.pql, D.star.inv){
  #generate m1 from t(0,D*)
  if(vars$m1>0) genData<-rmvt(ceiling(vars$m1/vars$no_cores),sigma=Dstarnotsparse,df=vars$zeta,type=c("shifted"))
  if(vars$m1==0) genData<-NULL		
  
  #generate m2 from N(u*,D*)
  if(m2>0) genData2<-genRand(vars$u.star,vars$D.star,ceiling(m2/vars$no_cores))
  if(m2==0) genData2<-NULL
  
  
  #generate m3 from N(u*,(Z'c''(Xbeta*+zu*)Z+D*^{-1})^-1)
  if(m3>0){
    Z=do.call(cbind,vars$mod.mcml$z)
    eta.star<-as.vector(vars$mod.mcml$x%*%beta.pql+Z%*%vars$u.star)
    if(vars$family.glmm$family.glmm=="bernoulli.glmm") {cdouble<-vars$family.glmm$cpp(eta.star)}
    if(vars$family.glmm$family.glmm=="poisson.glmm"){cdouble<-vars$family.glmm$cpp(eta.star)}
    if(vars$family.glmm$family.glmm=="binomial.glmm"){cdouble<-vars$family.glmm$cpp(eta.star, vars$ntrials)}
    #still a vector
    cdouble<-Diagonal(length(cdouble),cdouble)
    Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
    Sigmuh<-solve(Sigmuh.inv)
    genData3<-genRand(vars$u.star,Sigmuh,ceiling(m3/vars$no_cores))
  }
  if(m3==0) genData3<-NULL
  
  #	#these are from distribution based on data
  #	if(distrib=="tee")genData<-genRand(sigma.gen,s.pql,mod.mcml$z,m1,distrib="tee",gamm)
  #	if(distrib=="normal")genData<-genRand(sigma.pql,s.pql,mod.mcml$z,m1,distrib="normal",gamm)
  #	#these are from standard normal
  #	ones<-rep(1,length(sigma.pql))
  #	zeros<-rep(0,length(s.pql))
  #	genData2<-genRand(ones,zeros,mod.mcml$z,m2,distrib="normal",gamm)
  
  umat<-rbind(genData,genData2,genData3)
  m <- nrow(umat)
  list(umat=umat, m=m, Sigmuh.inv=Sigmuh.inv)
}

clusterSetRNGStream(vars$cl, 1234)

clusterExport(vars$cl, c("vars", "Dstarnotsparse", "m2", "m3", "beta.pql", "D.star.inv", "simulate", "genRand"), envir = environment())     #installing variables on each core
clusterEvalQ(vars$cl, umatparams <- simulate(vars=vars, Dstarnotsparse=Dstarnotsparse, m2=m2, m3=m3, beta.pql=beta.pql, D.star.inv=D.star.inv))

vars$nbeta <- 1

vars$p1=vars$p2=vars$p3=1/3

par<-c(6,1.5)
del<-rep(10^-8,2)

objfun<-glmm:::objfun

core2<-objfun(par=par, vars=vars)

umats <- clusterEvalQ(vars$cl, umatparams$umat)
umat <- Reduce(rbind, umats)

ms <- clusterEvalQ(vars$cl, umatparams$m)
m <- ms[[1]]

Sigmuh.invs <- clusterEvalQ(vars$cl, umatparams$Sigmuh.inv)
Sigmuh.inv <- Sigmuh.invs[[1]]

stopCluster(clust)

vars$cl <- makeCluster(1)

registerDoParallel(vars$cl)                   #making cluster usable with foreach
vars$no_cores <- length(vars$cl)

clusterExport(vars$cl, c("umat", "m", "Sigmuh.inv"), envir = environment())
clusterEvalQ(vars$cl, umatparams <- list(umat=umat, m=m, Sigmuh.inv=Sigmuh.inv))

core1<-objfun(par=par, vars=vars)

all.equal(core1, core2)

stopCluster(vars$cl)
