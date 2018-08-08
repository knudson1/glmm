#ntrials is a vector with length equal to length(y). if Bern or Poisson, ntrials is a vec of 1s

objfun <-
function(par, nbeta, nu.pql, umat, u.star, mod.mcml, family.glmm, cache, p1, p2, p3, m1, D.star, Sigmuh, Sigmuh.inv, zeta, ntrials, no_cores){

	beta<-par[1:nbeta]
	nu<-par[-(1:nbeta)]
	m<-nrow(umat)

	if (!missing(cache)) stopifnot(is.environment(cache))

	if(any(nu<=0)){
		out<-list(value=-Inf,gradient=rep(1,length(par)),hessian=as.matrix(c(rep(1,length(par)^2)),nrow=length(par)))
	return(out)
	}
	
	Z=do.call(cbind,mod.mcml$z)
	T<-length(mod.mcml$z)
	nrand<-lapply(mod.mcml$z,ncol)
	nrandom<-unlist(nrand)



	family.glmm<-getFamily(family.glmm)
	if(family.glmm$family.glmm=="bernoulli.glmm"){family_glmm=1}	
	if(family.glmm$family.glmm=="poisson.glmm"){family_glmm=2}	
	if(family.glmm$family.glmm=="binomial.glmm"){family_glmm=3}	

	Dstarinvdiag<-1/diag(D.star)
	D.star.inv<-diag(Dstarinvdiag)

	logdet.D.star.inv<-	-sum(log(diag(D.star)))
	logdet.Sigmuh.inv<-sum(log(eigen(Sigmuh.inv,symmetric=TRUE)$values))
 	myq<-nrow(D.star.inv)

	tconst<-tconstant(zeta,myq,Dstarinvdiag)

	#for the particular value of nu we're interested in, need to prep for distRandGenC
	eek<-getEk(mod.mcml$z)
	preDinvfornu<-Map("*",eek,(1/nu))
	Dinvfornu<-addVecs(preDinvfornu)
	logdetDinvfornu<-sum(log(Dinvfornu))
	Dinvfornu<-diag(Dinvfornu)
	
	meow<-rep(1,T+1)
	meow[1]<-0
	throwaway<-T+1
	meow[2:throwaway]<-cumsum(nrandom)
	
	pea<-c(p1,p2,p3)
	n<-nrow(mod.mcml$x)

##need to scale first m1 vectors of generated random effects by multiplying by A

#	preAfornu<-Map("*",eek,sqrt(nu))
#	Afornu<-addVecs(preAfornu)

#	for(k in 1:m1){
#		u.swoop<-umat[k,]
#		umat[k,]<-u.swoop*Afornu
#		}
	
	miniu <- NULL

	#parallelizing the calculations for the value of the log-likelihood approximation and gradient
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)
	clusterEvalQ(cl, library(itertools))
	clusterExport(cl, c("umat", "myq", "m", "mod.mcml", "n", "nbeta", "beta", "Z", "Dinvfornu", "logdetDinvfornu", "family_glmm", "D.star.inv", "logdet.D.star.inv", "u.star", "logdet.Sigmuh.inv", "pea", "T", "nrandom", "meow", "nu", "zeta", "tconst", "ntrials", "par", "Sigmuh.inv"), envir = environment())
	out <- foreach(miniu=isplitRows(umat, chunks = no_cores)) %dopar% {.C(C_valgrad, as.double(mod.mcml$y),as.double(t(miniu)), as.integer(myq), as.integer(nrow(miniu)), as.double(mod.mcml$x), as.integer(n), as.integer(nbeta), as.double(beta), as.double(Z), as.double(Dinvfornu), as.double(logdetDinvfornu),as.integer(family_glmm), as.double(D.star.inv), as.double(logdet.D.star.inv), as.double(u.star), as.double(Sigmuh.inv), as.double(logdet.Sigmuh.inv), pea=as.double(pea), nps=as.integer(length(pea)), T=as.integer(T), nrandom=as.integer(nrandom), meow=as.integer(meow),nu=as.double(nu), zeta=as.integer(zeta),tconst=as.double(tconst), v=double(m), ntrials=as.integer(ntrials), value=double(1),gradient=double(length(par)), b=double(nrow(miniu)))}
	stopCluster(cl)
	
	#combining the b's from each core into one vector
	b <- c()
	for(i in 1:no_cores){
	  b <- c(b, out[[i]]$b)
	}
	
	#recalculating the approximation of the log-likelihood value
	#finding a, the max of values
	vals <- c(rep(0, no_cores))
	for(i in 1:no_cores){
	  vals[i] <- out[[i]]$value
	}
	a <- max(vals)
	
	#storing normalized b[k], collected from values
	normalbk <- c(rep(0, no_cores))
	for(i in 1:no_cores){
	  normalbk[i] <- out[[i]][[4]]*exp(out[[i]]$value - a)
	}
	sumnbk <- sum(normalbk)
	
	value <- log(sumnbk/m) + a
	
	#recalculating the gradient
	gradient <- c(rep(0, length(out[[1]]$gradient)))
	for(j in 1: length(out[[1]]$gradient)){
	  gradadd <- 0
	  for(i in 1:length(normalbk)){
	    gradadd <- gradadd + normalbk[i]*out[[i]]$gradient[j]
	  }
	  gradient[j] <- gradadd/sumnbk
	}
	
	#parallelizing calculations for the hessian
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)
	clusterEvalQ(cl, library(itertools))
	clusterExport(cl, c("umat", "myq", "m", "mod.mcml", "n", "nbeta", "beta", "Z", "Dinvfornu", "logdetDinvfornu", "family_glmm", "D.star.inv", "logdet.D.star.inv", "u.star", "logdet.Sigmuh.inv", "pea", "T", "nrandom", "meow", "nu", "zeta", "tconst", "ntrials", "par", "Sigmuh.inv", "gradient"), envir = environment())
	out2 <- foreach(miniu=isplitRows(umat, chunks = no_cores)) %dopar% {.C(C_hess, as.double(mod.mcml$y),as.double(t(miniu)), as.integer(myq), as.integer(nrow(miniu)), as.double(mod.mcml$x), as.integer(n), as.integer(nbeta), as.double(beta), as.double(Z), as.double(Dinvfornu), as.double(logdetDinvfornu),as.integer(family_glmm), as.double(D.star.inv), as.double(logdet.D.star.inv), as.double(u.star), as.double(Sigmuh.inv), as.double(logdet.Sigmuh.inv), pea=as.double(pea), nps=as.integer(length(pea)), T=as.integer(T), nrandom=as.integer(nrandom), meow=as.integer(meow),nu=as.double(nu), zeta=as.integer(zeta),tconst=as.double(tconst), v=double(nrow(miniu)), ntrials=as.integer(ntrials),gradient=as.double(gradient),hessian=double((length(par))^2), b=as.double(b), length=as.integer(m))}
	stopCluster(cl)
	
	#adding hessian components
	hessian <- c(rep(0, length(out2[[1]]$hessian)))
	for(j in 1:length(out2[[1]]$hessian)){
	  hessadd <- 0
	  for(i in 1:no_cores){
	    hessadd <- hessadd + out2[[i]]$hessian[j]
	  }
	  hessian[j] <- hessadd
	}
	
	#making weights accessible
	weights <- c()
	for(i in 1:no_cores){
	  weights <- c(weights, out2[[i]]$v)
	}
  if (!missing(cache)) cache$weights<-weights		

  list(value=value,gradient=gradient,hessian=matrix(hessian,ncol=length(par),byrow=FALSE))
	
}

