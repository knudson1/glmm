#ntrials is a vector with length equal to length(y). if Bern or Poisson, ntrials is a vec of 1s

mcvcov <- function(object){

	mod <- object
	varcomps.equal <- mod$varcomps.equal
	levs<-ordered(unique(varcomps.equal))
	

	beta <- mod$beta
	nu <-mod$nu
	umat<-mod$umat
	m<-nrow(umat)

	x <- mod$x
	y <- mod$y
	random <- mod$z
	z<-list()
	for(i in 1:length(levs)){
	    if(levs[i]!=i) stop("The numbers in the vector varcomps.equal must be consecutive. You must start at 1 and then each entry must be the next consecutive number or a repeat of a previous number.")
	    these<-varcomps.equal==i
	    thesemats<-random[these]
	    z[[i]]<-do.call(cbind,thesemats) 
	}
	names(z)<-mod$varcomps.names #z is list of length = number of varcomps
	
	Z = do.call(cbind, z)
	Tee <-length(z)
	nrand<-lapply(mod$z, ncol)
	nrandom<-unlist(nrand)


	family.glmm<-getFamily(mod$family.glmm)
	if(family.glmm$family.glmm=="bernoulli.glmm"){family_glmm=1}	
	if(family.glmm$family.glmm=="poisson.glmm"){family_glmm=2}	
	if(family.glmm$family.glmm=="binomial.glmm"){family_glmm=3}	

	pvec <-mod$pvec
	p1 <- pvec[1]
	p2 <- pvec[2]
	p3 <- pvec[3]
	zeta <- mod$zeta
	beta.pql <- mod$beta.pql
	nu.pql <- mod$nu.pql
	
	#need to get D*
	eek<-getEk(mod$z)
	Aks<-Map("*",eek,nu.pql)
	D.star<-addVecs(Aks) 
	D.star<-diag(D.star)
	D.star.inv<-solve(D.star)

	#need to recreate the variance matrix of  imp sampling distribution
	
	eta.star<-as.vector(x%*%beta.pql+Z%*%mod$u.pql)
	cdouble<-bernoulli.glmm()$cpp(eta.star) #still a vector
	cdouble<-diag(cdouble)
	Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
	Sigmuh<-solve(Sigmuh.inv)

	Dstinvdiag<-diag(D.star.inv)

	logdet.D.star.inv<-	-sum(log(diag(D.star)))
	logdet.Sigmuh.inv<-sum(log(eigen(Sigmuh.inv,symmetric=TRUE)$values))
 	myq<-nrow(D.star.inv)

	tconst<-tconstant(zeta,myq,Dstinvdiag)

	#for the particular value of nu we're interested in, need to prep for distRandGenC
	eek<-getEk(mod$z)
	preDinvfornu<-Map("*",eek,(1/nu))
	Dinvfornu<-addVecs(preDinvfornu)
	logdetDinvfornu<-sum(log(Dinvfornu))
	Dinvfornu<-diag(Dinvfornu)
	
	meow<-rep(1,Tee+1)
	meow[1]<-0
	throwaway<-Tee+1
	meow[2:throwaway]<-cumsum(nrandom)
	
	pea<-mod$pvec
	n<-nrow(mod$x)

	nbeta<- length(beta)
	npar <- length(beta) + length(nu)
	squaretop <- rep(0,m)

	if(is.null(mod$weights)){
	  wts <- rep(1,length(y))
	} else{
	  wts <- mod$weights
	}

	# only YOU can prevent segfault errors
	# prevent segfault errors by checking dims
	stopifnot(length(y) == n)
	stopifnot(nrow(mod$x) == n)
	stopifnot(ncol(mod$x) == nbeta)
	stopifnot(length(beta) == nbeta)
	stopifnot(nrow(Z) == n)
	stopifnot(ncol(Z) == myq)
	stopifnot(nrow(Dinvfornu) == myq)
	stopifnot(ncol(Dinvfornu) == myq)
	stopifnot(length(logdetDinvfornu) == 1)
	stopifnot(length(family_glmm) == 1)
	stopifnot(nrow(D.star.inv) == myq)
	stopifnot(ncol(D.star.inv) == myq)
	stopifnot(length(logdet.D.star.inv) == 1)
    stopifnot(length(mod$u.pql) == myq)
    stopifnot(nrow(Sigmuh.inv) == myq)
    stopifnot(ncol(Sigmuh.inv) == myq)
    stopifnot(length(logdet.Sigmuh.inv) == 1)
    stopifnot(length(Tee) == 1)
    stopifnot(length(nrandom) == Tee)
    stopifnot(length(nrandom) == length(nu))
    stopifnot(length(meow) == (Tee+1))
    stopifnot(length(nu) == Tee)
    stopifnot(length(zeta) == 1)
    stopifnot(length(tconst) == 1)
    stopifnot(length(mod$mod.mcml$ntrials) == n)
    stopifnot(length(wts) == n)
    
	
	stuff<-.C(C_mcsec, 
	          as.double(0.0), # gamma: scalar
	          as.double(0.0), # thing: scalar
	          as.double(squaretop), # squaretop: vector length m
	          numsum = as.double(rep(0, npar^2)), # numsum: vector of length  npar^2
	          as.double(mod$y), # y: vector of length n
	          as.double(t(umat)), # Umat: myq by m matrix.
	          as.integer(myq), # scalar.
	          as.integer(m), # scalar.
	          as.double(mod$x), # n by nbeta matrix
	          as.integer(n), #scalar
	          as.integer(nbeta),  #scalar
	          as.double(beta), # vector length nbeta
	          as.double(Z), # n by myq matrix
	          as.double(Dinvfornu), # myq x myq
	          as.double(logdetDinvfornu), #scalar
	          as.integer(family_glmm), #scalar
	          as.double(D.star.inv), # myq x myq.
	          as.double(logdet.D.star.inv),  #scalar
	          as.double(mod$u.pql), # length = myq
	          as.double(Sigmuh.inv),  # myq x myq.
	          as.double(logdet.Sigmuh.inv), # scalar
	          pea=as.double(pea), # vector of length nps
	          nps=as.integer(length(pea)), #scalar
	          Tee=as.integer(Tee), # scalar equal to the number of variance components
	          nrandom=as.integer(nrandom), # vector of length T. 
	          meow=as.integer(meow), # vector of length T+1.
	          nu=as.double(nu), # vector of length T
	          zeta=as.integer(zeta), #scalar
	          tconst=as.double(tconst), #scalar.
	          ntrials=as.integer(mod$mod.mcml$ntrials), #vector of length n
	          as.double(0.0), # lfuval: scalar
	          as.double(0.0), #  lfyuval: scalar
	          wts=as.double(wts)) # vector length n. (to calculate weighted likelihood)

	vhatnum <- (1/m)*stuff[[4]]
	vhatdenom <- ( stuff[[1]]   )^2

	vhatvec <- vhatnum/vhatdenom
	Vhat <- matrix(vhatvec, nrow=npar)

	Uhat <- mod$loglike.hessian
	Uhatinv <- qr.solve(Uhat)

	out <- Uhatinv %*% Vhat %*% Uhatinv/m

	cf<-c(coef(object),varcomps(object))
	pnames<-names(cf)
	rownames(out) <- pnames
	colnames(out) <- pnames

	out

}

#ntrials is a vector with length equal to length(y). if Bern or Poisson, ntrials is a vec of 1s

mcse <- function(object){


	UinvVUinv <- mcvcov(object)

	MCSE <- sqrt(diag(UinvVUinv))
	MCSE

}




