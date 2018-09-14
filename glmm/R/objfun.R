#ntrials is a vector with length equal to length(y). if Bern or Poisson, ntrials is a vec of 1s

objfun <-
function(par, cache, vars){
  
  vars$par <- par
  
	vars$beta<-beta<-vars$par[1:vars$nbeta]
	vars$nu<-vars$par[-(1:vars$nbeta)]
	vars$m<-nrow(vars$umat)

	if (!missing(cache)) stopifnot(is.environment(cache))

	if(any(vars$nu<=0)){
		out<-list(value=-Inf,gradient=rep(1,length(vars$par)),hessian=as.matrix(c(rep(1,length(vars$par)^2)),nrow=length(vars$par)))
	return(out)
	}
	
	vars$Z=do.call(cbind,vars$mod.mcml$z)
	vars$T<-length(vars$mod.mcml$z)
	nrand<-lapply(vars$mod.mcml$z,ncol)
	vars$nrandom<-unlist(nrand)



	vars$family.glmm<-getFamily(vars$family.glmm)
	if(vars$family.glmm$family.glmm=="bernoulli.glmm"){vars$family_glmm=1}	
	if(vars$family.glmm$family.glmm=="poisson.glmm"){vars$family_glmm=2}	
	if(vars$family.glmm$family.glmm=="binomial.glmm"){vars$family_glmm=3}	

	Dstarinvdiag<-1/diag(vars$D.star)
	vars$D.star.inv<-diag(Dstarinvdiag)

	vars$logdet.D.star.inv<-	-sum(log(diag(vars$D.star)))
	vars$logdet.Sigmuh.inv<-sum(log(eigen(vars$Sigmuh.inv,symmetric=TRUE)$values))
 	vars$myq<-nrow(vars$D.star.inv)

	vars$tconst<-tconstant(vars$zeta,vars$myq,Dstarinvdiag)

	#for the particular value of nu we're interested in, need to prep for distRandGenC
	eek<-getEk(vars$mod.mcml$z)
	preDinvfornu<-Map("*",eek,(1/vars$nu))
	vars$Dinvfornu<-addVecs(preDinvfornu)
	vars$logdetDinvfornu<-sum(log(vars$Dinvfornu))
	vars$Dinvfornu<-diag(vars$Dinvfornu)
	
	vars$meow<-rep(1,vars$T+1)
	vars$meow[1]<-0
	throwaway<-vars$T+1
	vars$meow[2:throwaway]<-cumsum(vars$nrandom)
	
	vars$pea<-c(vars$p1,vars$p2,vars$p3)
	vars$n<-nrow(vars$mod.mcml$x)

##need to scale first m1 vectors of generated random effects by multiplying by A

#	preAfornu<-Map("*",eek,sqrt(nu))
#	Afornu<-addVecs(preAfornu)

#	for(k in 1:m1){
#		u.swoop<-umat[k,]
#		umat[k,]<-u.swoop*Afornu
#		}
	
	miniu <- NULL
	minib <- NULL

	#parallelizing the calculations for the value of the log-likelihood approximation and gradient
	cl <- makeCluster(vars$no_cores)              #making cluster environment
	registerDoParallel(cl)                   #making cluster usable with foreach
	clusterEvalQ(cl, library(itertools))     #installing itertools library on each core
	clusterExport(cl, c("umat", "myq", "m", "mod.mcml", "n", "nbeta", "beta", "Z", "Dinvfornu", "logdetDinvfornu", "family_glmm", "D.star.inv", "logdet.D.star.inv", "u.star", "logdet.Sigmuh.inv", "pea", "T", "nrandom", "meow", "nu", "zeta", "tconst", "ntrials", "par", "Sigmuh.inv"), envir = vars)     #installing variables on each core
	out <- foreach(miniu=isplitRows(vars$umat, chunks = vars$no_cores)) %dopar% {.C(C_valgrad, as.double(mod.mcml$y),as.double(t(miniu)), as.integer(myq), as.integer(nrow(miniu)), as.double(mod.mcml$x), as.integer(n), as.integer(nbeta), as.double(beta), as.double(Z), as.double(Dinvfornu), as.double(logdetDinvfornu),as.integer(family_glmm), as.double(D.star.inv), as.double(logdet.D.star.inv), as.double(u.star), as.double(Sigmuh.inv), as.double(logdet.Sigmuh.inv), pea=as.double(pea), nps=as.integer(length(pea)), T=as.integer(T), nrandom=as.integer(nrandom), meow=as.integer(meow),nu=as.double(nu), zeta=as.integer(zeta),tconst=as.double(tconst), v=double(m), ntrials=as.integer(ntrials), value=double(1),gradient=double(length(par)), b=double(nrow(miniu)))}
	#foreach runs loop in parallel, dopar operator sends each chunk of umat to seperate core and runs the .C function
	stopCluster(cl)     #closing the cluster environment
	
	#combining the b's from each core into one vector
	bs <- as.list(rep(0, vars$no_cores))
	for(i in 1:vars$no_cores){
	  bs[[i]] <- out[[i]]$b
	}
	
	vars$b <- Reduce(c, bs)
	
	#recalculating the approximation of the log-likelihood value
	#finding a, the max of values
	vals <- c(rep(0, vars$no_cores))
	for(i in 1:vars$no_cores){
	  vals[i] <- out[[i]]$value
	}
	a <- max(vals)
	
	#storing normalized b[k], collected from values
	normalbk <- c(rep(0, vars$no_cores))
	for(i in 1:vars$no_cores){
	  normalbk[i] <- out[[i]][[4]]*exp(out[[i]]$value - a)
	}
	sumnbk <- sum(normalbk)
	
	value <- log(sumnbk/vars$m) + a
	
	#recalculating the gradient
	vars$gradient <- c(rep(0, length(out[[1]]$gradient)))
	for(j in 1: length(out[[1]]$gradient)){
	  gradadd <- 0
	  for(i in 1:length(normalbk)){
	    gradadd <- gradadd + normalbk[i]*out[[i]]$gradient[j]
	  }
	  vars$gradient[j] <- gradadd/sumnbk
	}
	
	#parallelizing calculations for the hessian
	cl <- makeCluster(vars$no_cores)
	registerDoParallel(cl)
	clusterEvalQ(cl, library(itertools))
	clusterExport(cl, c("umat", "myq", "m", "mod.mcml", "n", "nbeta", "beta", "Z", "Dinvfornu", "logdetDinvfornu", "family_glmm", "D.star.inv", "logdet.D.star.inv", "u.star", "logdet.Sigmuh.inv", "pea", "T", "nrandom", "meow", "nu", "zeta", "tconst", "ntrials", "par", "Sigmuh.inv", "gradient", "b"), envir = vars)
	out2 <- foreach(miniu=isplitRows(vars$umat, chunks = vars$no_cores), minib=isplitVector(vars$b, chunks = vars$no_cores)) %dopar% {.C(C_hess, as.double(mod.mcml$y),as.double(t(miniu)), as.integer(myq), as.integer(nrow(miniu)), as.double(mod.mcml$x), as.integer(n), as.integer(nbeta), as.double(beta), as.double(Z), as.double(Dinvfornu), as.double(logdetDinvfornu),as.integer(family_glmm), as.double(D.star.inv), as.double(logdet.D.star.inv), as.double(u.star), as.double(Sigmuh.inv), as.double(logdet.Sigmuh.inv), pea=as.double(pea), nps=as.integer(length(pea)), T=as.integer(T), nrandom=as.integer(nrandom), meow=as.integer(meow),nu=as.double(nu), zeta=as.integer(zeta),tconst=as.double(tconst), v=double(nrow(miniu)), ntrials=as.integer(ntrials),gradient=as.double(gradient),hessian=double((length(par))^2), b=as.double(b), length=as.integer(m), q=as.double(minib))}
	stopCluster(cl)
	
	#adding hessian components
	hessian <- c(rep(0, length(out2[[1]]$hessian)))
	for(j in 1:length(out2[[1]]$hessian)){
	  hessadd <- 0
	  for(i in 1:vars$no_cores){
	    hessadd <- hessadd + out2[[i]]$hessian[j]
	  }
	  hessian[j] <- hessadd
	}
	
	#making weights accessible
	weight <- as.list(rep(0, vars$no_cores))
	for(i in 1:vars$no_cores){
	  weight[[i]] <- out2[[i]]$v
	}
	weights <- Reduce(c, weight)
  if (!missing(cache)) cache$weights<-weights		

  list(value=value,gradient=vars$gradient,hessian=matrix(hessian,ncol=length(vars$par),byrow=FALSE))
	
}

