library(glmm)
data(BoothHobert)

clust <- makeCluster(1)
set.seed(1234)
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), data=BoothHobert, family.glmm=bernoulli.glmm, m=21, doPQL=TRUE, debug=TRUE, cluster=clust)

vars1 <- new.env(parent = emptyenv())

mod.mcml<-mod.mcml1$mod.mcml
#z<-mod.mcml$z[[1]]
#x<-mod.mcml$x
#y<-mod.mcml$y

y<- mod.mcml$y[1:149]
x<-mod.mcml$x[1:149]
z<- mod.mcml$z[[1]][1:149]

if(is.null(mod.mcml1$weights)){
  wts <- rep(1, length(y))
} else{
  wts <- mod.mcml1$weights
}

stuff<-mod.mcml1$debug
beta.pql<-stuff$beta.pql
nu.pql<-stuff$nu.pql
vars1$nu.pql<-nu.pql
u.pql<-u.star<-stuff$u.star
vars1$u.star<-u.star
umat<-stuff$umat
vars1$family.glmm<-mod.mcml1$family.glmm
vars1$umat<-stuff$umat
vars1$newm <- nrow(vars1$umat)
vars1$ntrials<-1
D.star.inv <- Dstarnotsparse <- vars1$D.star <- as.matrix(stuff$D.star)

vars1$m1 <- stuff$m1
m2 <- stuff$m2
m3 <- stuff$m3

vars1$zeta <- 5

family.glmm<-bernoulli.glmm

objfun<-glmm:::objfun
getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

vars1$cl <- mod.mcml1$cluster
registerDoParallel(vars1$cl)                   #making cluster usable with foreach
vars1$no_cores <- length(vars1$cl)

vars1$mod.mcml<-mod.mcml1$mod.mcml

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

clusterSetRNGStream(vars1$cl, 1234)

clusterExport(vars1$cl, c("vars1", "Dstarnotsparse", "m2", "m3", "beta.pql", "D.star.inv", "simulate", "genRand"), envir = environment())     #installing variables on each core
noprint <- clusterEvalQ(vars1$cl, umatparams <- simulate(vars=vars1, Dstarnotsparse=Dstarnotsparse, m2=m2, m3=m3, beta.pql=beta.pql, D.star.inv=D.star.inv))

vars1$nbeta <- 1

vars1$p1=vars1$p2=vars1$p3=1/3

vars1$wts <- c(rep(1,149),2)

par<-c(6,1.5)
del<-rep(10^-9,2)

objfun<-glmm:::objfun

objfun1<-objfun(par=par, vars=vars1)

umats <- clusterEvalQ(vars1$cl, umatparams$umat)
umatT <- Reduce(rbind, umats)  #umatT is the total umat

Sigmuh.invs <- clusterEvalQ(vars1$cl, umatparams$Sigmuh.inv)
Sigmuh.invT <- Sigmuh.invs[[1]]

stopCluster(clust)

############################################
#this should be the same as el
getFamily<-glmm:::getFamily
elR <-
  function(Y,X,eta,family.mcml,wts){
    family.mcml<-getFamily(family.mcml)
    neta<-length(eta)
    
    if(family.mcml$family.glmm=="bernoulli.glmm"){
      foo<-.C(glmm:::C_cum3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(1), ntrials=as.integer(1), wts=as.double(wts), cumout=double(1))$cumout
      mu<-.C(glmm:::C_cp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(1), ntrials=as.integer(1), cpout=double(neta))$cpout
      cdub<-.C(glmm:::C_cpp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(1), ntrials=as.integer(1), cppout=double(neta))$cppout
    }
    if(family.mcml$family.glmm=="poisson.glmm"){
      foo<-.C(glmm:::C_cum3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(2), ntrials=as.integer(1), wts=as.double(wts), cumout=double(1))$cumout
      mu<-.C(glmm:::C_cp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),ntrials=as.integer(1),cpout=double(neta))$cpout
      cdub<-.C(glmm:::C_cpp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),ntrials=as.integer(1),cppout=double(neta))$cppout
    }
    
    value<-as.numeric(Y%*%eta-foo)
    gradient<-t(X)%*%(Y-mu)	
    cdubmat<-diag(cdub)
    hessian<-t(X)%*%(-cdubmat)%*%X
    
    list(value=value,gradient=gradient,hessian=hessian)
  }

NEWelR <-
  function(Y,X,eta,family.mcml,wts){
    family.mcml<-getFamily(family.mcml)
    neta<-length(eta)
    
    if(family.mcml$family.glmm=="bernoulli.glmm"){
      foo<-.C(glmm:::C_cum3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(1), ntrials=as.integer(1), wts=as.double(wts), cumout=double(1))$cumout
      mu<-.C(glmm:::C_cp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(1), ntrials=as.integer(1), cpout=double(neta))$cpout
      cdub<-.C(glmm:::C_cpp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(1), ntrials=as.integer(1), cppout=double(neta))$cppout
    }
    if(family.mcml$family.glmm=="poisson.glmm"){
      foo<-.C(glmm:::C_cum3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(2), ntrials=as.integer(1), wts=as.double(wts), cumout=double(1))$cumout
      mu<-.C(glmm:::C_cp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),ntrials=as.integer(1),cpout=double(neta))$cpout
      cdub<-.C(glmm:::C_cpp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),ntrials=as.integer(1),cppout=double(neta))$cppout
    }
    
    wtsmat <- diag(wts)
    wtX <- wtsmat%*%X
    
    value<-as.numeric(Y%*%wtsmat%*%eta-foo)
    gradient<-t(wtX)%*%(Y-mu)	
    cdubmat<-diag(cdub)
    hessian<-t(wtX)%*%(-cdubmat)%*%X
    
    list(value=value,gradient=gradient,hessian=hessian)
  }

#compare elR and NEWelR for weights all equal 1
eta1<-rep(2,150)
thatALL1<-elR(mod.mcml$y,mod.mcml$x,eta1,family.mcml=bernoulli.glmm, wts=wts) 
thisALL1 <- NEWelR(mod.mcml$y,mod.mcml$x,eta1,family.mcml=bernoulli.glmm, wts=wts)
all.equal(as.numeric(thatALL1$value),as.numeric(thisALL1$value))
all.equal(as.numeric(thatALL1$gradient),as.numeric(thisALL1$gradient))
all.equal(as.numeric(thatALL1$hessian),as.numeric(thisALL1$hessian))

#compare NEWelR and elc for weights all equal 1
thoseALL1<-.C(glmm:::C_elc, as.double(mod.mcml$y), as.double(mod.mcml$x), as.integer(nrow(mod.mcml$x)), as.integer(ncol(mod.mcml$x)), as.double(eta1), as.integer(1), as.integer(1), wts=as.double(rep(1,150)), value=double(1), gradient=double(ncol(mod.mcml$x)), hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.numeric(thoseALL1$value),as.numeric(thisALL1$value))
all.equal(as.numeric(thoseALL1$gradient),as.numeric(thisALL1$gradient))
all.equal(as.numeric(thoseALL1$hessian),as.numeric(thisALL1$hessian))

#finite differences for NEWelR
del<- 10^-9
thisdel <- NEWelR(mod.mcml$y,mod.mcml$x,eta1+del,family.mcml=bernoulli.glmm, wts=wts) 

all.equal(as.vector(thisALL1$gradient*del),thisdel$value-thisALL1$value)
all.equal(as.vector(thisALL1$hessian*del),as.vector(thisdel$gradient-thisALL1$gradient))

#finite differences for elc
thosedel <- .C(glmm:::C_elc, as.double(mod.mcml$y), as.double(mod.mcml$x), as.integer(nrow(mod.mcml$x)), as.integer(ncol(mod.mcml$x)), as.double(eta1+del), as.integer(1), as.integer(1), wts=as.double(rep(1,150)), value=double(1), gradient=double(ncol(mod.mcml$x)), hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.vector(thoseALL1$gradient*del),thosedel$value-thoseALL1$value)
all.equal(as.vector(thoseALL1$hessian*del),as.vector(thosedel$gradient-thoseALL1$gradient))

#compare elc to elval
elvalout<-.C(glmm:::C_elval, as.double(mod.mcml$y), as.integer(nrow(mod.mcml$x)), as.integer(ncol(mod.mcml$x)), as.double(eta1), as.integer(1), as.integer(1), wts=as.double(rep(1,150)), value=double(1))
all.equal(as.numeric(thoseALL1$value),elvalout$value)

#compare elc to elGH
elGHout<-.C(glmm:::C_elGH,as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta1),as.integer(1), as.integer(1), wts=as.double(rep(1,150)), gradient=double(ncol(mod.mcml$x)),hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.numeric(thoseALL1$gradient),elGHout$gradient)
all.equal(as.numeric(thoseALL1$hessian),elGHout$hessian)

#compare elR with duplicate last entry and NEWelR with wts=2 for last entry
BoothHobertDub <- rbind(BoothHobert, BoothHobert[nrow(BoothHobert),])

clust <- makeCluster(1)
eta2<-rep(2,151)
set.seed(1234)
mod.mcml2<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), data=BoothHobertDub, family.glmm=bernoulli.glmm, m=21, doPQL=TRUE, debug=TRUE, cluster=clust)

vars2 <- new.env(parent = emptyenv())

#z<-mod.mcml2$z[[1]]
#x<-mod.mcml2$x
#y<-mod.mcml2$y

y<-mod.mcml2$y[1:149]
x <- mod.mcml2$x[1:149]
z<- mod.mcml2$z[[1]][1:149,]

mod.mcml<-list(x=x, y=y, z=list(z), ntrials=rep(1,149))

if(is.null(mod.mcml2$weights)){
  wts <- rep(1, length(y))
} else{
  wts <- mod.mcml2$weights
}

stuff<-mod.mcml2$debug
beta.pql<-stuff$beta.pql
nu.pql<-stuff$nu.pql
vars2$nu.pql<-nu.pql
u.pql<-u.star<-stuff$u.star
vars2$u.star<-u.star
umat<-stuff$umat
vars2$family.glmm<-mod.mcml2$family.glmm
vars2$umat<-stuff$umat
vars2$newm <- nrow(vars2$umat)
vars2$ntrials<-1
D.star.inv <- Dstarnotsparse <- vars2$D.star <- as.matrix(stuff$D.star)

vars2$m1 <- stuff$m1
m2 <- stuff$m2
m3 <- stuff$m3

vars2$zeta <- 5

vars2$cl <- mod.mcml2$cluster
registerDoParallel(vars2$cl)                   #making cluster usable with foreach
vars2$no_cores <- length(vars2$cl)

#vars2$mod.mcml<-mod.mcml2$mod.mcml
vars2$mod.mcml<-mod.mcml

getEk<-glmm:::getEk
addVecs<-glmm:::addVecs
genRand<-glmm:::genRand

clusterExport(vars2$cl, c("umatT", "vars2", "Sigmuh.invT"), envir = environment())
noprint <- clusterEvalQ(vars2$cl, umatparams <- list(umat=umatT, m=nrow(umatT), Sigmuh.inv=Sigmuh.invT))

vars2$nbeta <- 1

vars2$p1=vars2$p2=vars2$p3=1/3

vars2$wts <- as.vector(wts)

objfun2<-objfun(par=par, vars=vars2)

stopCluster(clust)

this2<-NEWelR(mod.mcml$y,mod.mcml$x,eta1,family.mcml=bernoulli.glmm, wts=c(rep(1,149),2)) 
that2 <- elR(mod.mcml2$y,mod.mcml2$x,eta2,family.mcml=bernoulli.glmm, wts=wts)
those2 <- .C(glmm:::C_elc, as.double(mod.mcml$y), as.double(mod.mcml$x), as.integer(nrow(mod.mcml$x)), as.integer(ncol(mod.mcml$x)), as.double(eta1), as.integer(1), as.integer(1), wts=as.double(c(rep(1,149),2)), value=double(1), gradient=double(ncol(mod.mcml$x)), hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.numeric(that2$value),as.numeric(this2$value))
all.equal(as.numeric(that2$gradient),as.numeric(this2$gradient))
all.equal(as.numeric(that2$hessian),as.numeric(this2$hessian))

#comparing NEWelR C version
all.equal(as.numeric(those2$value),as.numeric(this2$value))
all.equal(as.numeric(those2$gradient),as.numeric(this2$gradient))
all.equal(as.numeric(those2$hessian),as.numeric(this2$hessian))

#comparing objfun output
all.equal(objfun1,objfun2)
