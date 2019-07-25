library(glmm)
data(BoothHobert)

set.seed(1234)
#model with all weights at 1, no duplicate data points in data set
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), data=BoothHobert, family.glmm=bernoulli.glmm, m=10^2, doPQL=TRUE, debug=TRUE)

#weights are determined from model (should be all 1)
if(is.null(mod.mcml1$weights)){
  wts <- rep(1, length(mod.mcml1$y))
} else{
  wts <- mod.mcml1$weights
}

############################################
getFamily<-glmm:::getFamily
#el without weights (in R)
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

#el with weights (in R)
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

########################################################
#compare elR and NEWelR for weights all equal 1: to make sure elR and NEWelR work the same with no weighting scheme
eta1<-rep(2,150)
mod.mcml<-mod.mcml1
thatALL1<-elR(mod.mcml$y,mod.mcml$x,eta1,family.mcml=bernoulli.glmm, wts=wts) 
thisALL1 <- NEWelR(mod.mcml$y,mod.mcml$x,eta1,family.mcml=bernoulli.glmm, wts=wts)
all.equal(as.numeric(thatALL1$value),as.numeric(thisALL1$value))
all.equal(as.numeric(thatALL1$gradient),as.numeric(thisALL1$gradient))
all.equal(as.numeric(thatALL1$hessian),as.numeric(thisALL1$hessian))

#compare NEWelR and elc for weights all equal 1: to make sure elc and NEWelR work the same with no weighting scheme
thoseALL1<-.C(glmm:::C_elc, as.double(mod.mcml$y), as.double(mod.mcml$x), as.integer(nrow(mod.mcml$x)), as.integer(ncol(mod.mcml$x)), as.double(eta1), as.integer(1), as.integer(1), wts=as.double(rep(1,150)), value=double(1), gradient=double(ncol(mod.mcml$x)), hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.numeric(thoseALL1$value),as.numeric(thisALL1$value))
all.equal(as.numeric(thoseALL1$gradient),as.numeric(thisALL1$gradient))
all.equal(as.numeric(thoseALL1$hessian),as.numeric(thisALL1$hessian))

#finite differences for NEWelR, weights all 1
del<- 10^-9
thisdel <- NEWelR(mod.mcml$y,mod.mcml$x,eta1+del,family.mcml=bernoulli.glmm, wts=wts) 

all.equal(as.vector(thisALL1$gradient*del),thisdel$value-thisALL1$value)
all.equal(as.vector(thisALL1$hessian*del),as.vector(thisdel$gradient-thisALL1$gradient))

#finite differences for elc, weights all 1
thosedel <- .C(glmm:::C_elc, as.double(mod.mcml$y), as.double(mod.mcml$x), as.integer(nrow(mod.mcml$x)), as.integer(ncol(mod.mcml$x)), as.double(eta1+del), as.integer(1), as.integer(1), wts=as.double(rep(1,150)), value=double(1), gradient=double(ncol(mod.mcml$x)), hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.vector(thoseALL1$gradient*del),thosedel$value-thoseALL1$value)
all.equal(as.vector(thoseALL1$hessian*del),as.vector(thosedel$gradient-thoseALL1$gradient))

#compare elc to elval, weights all 1: value should be the same
elvalout<-.C(glmm:::C_elval, as.double(mod.mcml$y), as.integer(nrow(mod.mcml$x)), as.integer(ncol(mod.mcml$x)), as.double(eta1), as.integer(1), as.integer(1), wts=as.double(rep(1,150)), value=double(1))
all.equal(as.numeric(thoseALL1$value),elvalout$value)

#compare elc to elGH, weights all 1: gradient and hessian should be the same
elGHout<-.C(glmm:::C_elGH,as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta1),as.integer(1), as.integer(1), wts=as.double(rep(1,150)), gradient=double(ncol(mod.mcml$x)),hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.numeric(thoseALL1$gradient),elGHout$gradient)
all.equal(as.numeric(thoseALL1$hessian),elGHout$hessian)

#BoothHobert with 151 data points instead of 150 (150th data point duplicated)
BoothHobertDub <- rbind(BoothHobert, BoothHobert[nrow(BoothHobert),])

eta2<-rep(2,151)
set.seed(1234)
#model using duplicated data, all weights are 1
mod.mcml2<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), data=BoothHobertDub, family.glmm=bernoulli.glmm, m=10^2, doPQL=TRUE, debug=TRUE)

#151 weights (all 1)
if(is.null(mod.mcml2$weights)){
  wts <- rep(1, length(mod.mcml2$y))
} else{
  wts <- mod.mcml2$weights
}

#compare elR with BoothHobertDub and all weights 1 versus NEWelR with BoothHobert and first 149 wights 1 and weight 150 as 2
this2<-NEWelR(mod.mcml$y,mod.mcml$x,eta1,family.mcml=bernoulli.glmm, wts=c(rep(1,149),2)) 
that2 <- elR(mod.mcml2$y,mod.mcml2$x,eta2,family.mcml=bernoulli.glmm, wts=wts)
all.equal(as.numeric(that2$value),as.numeric(this2$value))
all.equal(as.numeric(that2$gradient),as.numeric(this2$gradient))
all.equal(as.numeric(that2$hessian),as.numeric(this2$hessian))

#compare NEWelR with BoothHobert and first 149 wights 1 and weight 150 as 2 versus elc with BoothHobert and first 149 wights 1 and weight 150 as 2
those2 <- .C(glmm:::C_elc, as.double(mod.mcml$y), as.double(mod.mcml$x), as.integer(nrow(mod.mcml$x)), as.integer(ncol(mod.mcml$x)), as.double(eta1), as.integer(1), as.integer(1), wts=as.double(c(rep(1,149),2)), value=double(1), gradient=double(ncol(mod.mcml$x)), hessian=double((ncol(mod.mcml$x)^2)))
all.equal(as.numeric(those2$value),as.numeric(this2$value))
all.equal(as.numeric(those2$gradient),as.numeric(this2$gradient))
all.equal(as.numeric(those2$hessian),as.numeric(this2$hessian))
