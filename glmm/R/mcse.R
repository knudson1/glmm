#ntrials is a vector with length equal to length(y). if Bern or Poisson, ntrials is a vec of 1s

mcse <- function(mod){

	getFamily <- glmm:::getFamily
	getEk<-glmm:::getEk
	addVecs<-glmm:::addVecs
	tconstant <- glmm:::tconstant
	addMats<-glmm:::addMats
	

	beta <- mod$beta
	nu <-mod$nu
	umat<-mod$umat
	m<-nrow(umat)

	x <- mod$x
	y <- mod$y
	Z=do.call(cbind,mod$z)
	T<-length(mod$z)
	nrand<-lapply(mod$z,ncol)
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
	
	meow<-rep(1,T+1)
	meow[1]<-0
	throwaway<-T+1
	meow[2:throwaway]<-cumsum(nrandom)
	
	pea<-mod$pvec
	n<-nrow(mod$x)

	nbeta<- length(beta)

	squaretop <- rep(0,m)

	

	stuff<-.C("mcsec", as.double(0.0), as.double(0.0), as.double(squaretop), as.double(mod$y),as.double(t(umat)), as.integer(myq), as.integer(m), as.double(mod$x), as.integer(n), as.integer(nbeta), as.double(beta), as.double(Z), as.double(Dinvfornu), as.double(logdetDinvfornu),as.integer(family_glmm), as.double(D.star.inv), as.double(logdet.D.star.inv), as.double(mod$u.pql), as.double(Sigmuh.inv), as.double(logdet.Sigmuh.inv), pea=as.double(pea), nps=as.integer(length(pea)), T=as.integer(T), nrandom=as.integer(nrandom), meow=as.integer(meow),nu=as.double(nu), zeta=as.integer(zeta),tconst=as.double(tconst), v=double(m), ntrials=as.integer(mod$mod.mcml$ntrials), value=double(1),gradient=double(length(par)),hessian=double((length(par))^2),PACKAGE="glmm")






}

