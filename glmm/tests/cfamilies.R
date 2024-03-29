
library(glmm)

#First, define checks for bernoulli
cum <- function(eta){ 
	if(eta>0) eta+log1p(exp(-eta))
	else log1p(exp(eta))}
cum<-Vectorize(cum)
cp <- function(eta) {1/(1+exp(-eta))}
cpp<-function(eta) {(1/(1+exp(-eta)))*(1/(1+exp(eta)))}


#check cum3 for a vector containing both positive and negative values of eta
eta<-seq(-1,1,.1)
neta<-length(eta)
ntrials <- wts<-rep(1,neta)
that<-.C(glmm:::C_cum3,eta=as.double(eta),neta=as.integer(neta), type=as.integer(1), ntrials=as.integer(ntrials), wts=as.double(wts),cumout=double(1))
all.equal(sum(cum(eta)),that$cumout)

#check cp3
that<-.C(glmm:::C_cp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(1), ntrials=as.integer(ntrials), cpout=double(neta))
all.equal(cp(eta),that$cpout)

#check cpp3
that<-.C(glmm:::C_cpp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(1), ntrials=as.integer(ntrials), cppout=double(neta))
all.equal(cpp(eta),that$cppout)

######################
#now move onto Poisson
cum <- function(eta) exp(eta)
cp <- function(eta) exp(eta)
cpp<-function(eta) exp(eta)
eta<-c(4.5,5,5.5)
neta<-length(eta)
ntrials <- wts<-rep(1,neta)

#check cum3
that<-.C(glmm:::C_cum3, eta=as.double(eta), neta=as.integer(neta), type=as.integer(2), ntrials=as.integer(ntrials), wts=as.double(wts), cumout=double(1))
all.equal(sum(cum(eta)),that$cumout)

#check cp3
that<-.C(glmm:::C_cp3,eta=as.double(eta), neta=as.integer(neta), type=as.integer(2), ntrials=as.integer(ntrials), cpout=double(neta))
all.equal(cp(eta), that$cpout)

#check cpp3
that<-.C(glmm:::C_cpp3,eta=as.double(eta),neta=as.integer(neta),type=as.integer(2), ntrials=as.integer(ntrials), cppout=double(neta))
all.equal(cpp(eta), that$cppout)

