library(glmm)
data(BoothHobert)

set.seed(1234)
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), data=BoothHobert, family.glmm=bernoulli.glmm, m=21, doPQL=TRUE, debug=TRUE)

mod.mcml<-mod.mcml1$mod.mcml
z<-mod.mcml$z[[1]]
x<-mod.mcml$x
y<-mod.mcml$y

if(is.null(mod.mcml1$weights)){
  wts <- rep(1, length(y))
} else{
  wts <- mod.mcml1$weights
}

stuff<-mod.mcml1$debug
beta.pql<-stuff$beta.pql
nu.pql<-stuff$nu.pql
u.pql<-u.star<-stuff$u.star
umat<-stuff$umat

family.glmm<-bernoulli.glmm

objfun<-glmm:::objfun
getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

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

eta2<-rep(2,151)
set.seed(1234)
mod.mcml2<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), data=BoothHobertDub, family.glmm=bernoulli.glmm, m=21, doPQL=TRUE, debug=TRUE)

z<-mod.mcml2$z[[1]]
x<-mod.mcml2$x
y<-mod.mcml2$y

if(is.null(mod.mcml2$weights)){
  wts <- rep(1, length(mod.mcml2$y))
} else{
  wts <- mod.mcml2$weights
}

stuff<-mod.mcml2$debug
beta.pql<-stuff$beta.pql
nu.pql<-stuff$nu.pql
u.pql<-u.star<-stuff$u.star
umat<-stuff$umat

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