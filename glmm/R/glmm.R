utils::globalVariables("umatparams")

glmm <-
  function(fixed,random, varcomps.names,data, family.glmm, m,varcomps.equal, weights = NULL, doPQL=TRUE, debug=FALSE,p1=1/3,p2=1/3,p3=1/3,rmax=1000,iterlim=1000,par.init=NULL,zeta=5, cluster=NULL)
  {
    if(missing(varcomps.names)) stop("Names for the variance components must be supplied through varcomps.names")
    if(is.vector(varcomps.names)!=1) stop("varcomps.names must be a vector")
    if(missing(varcomps.equal)){
      varcomps.equal<- c(1:length(varcomps.names))}
    call<-match.call()
    
    #vars will store all variables needed for valgrad and hess in objfun, and trust
    vars <- new.env(parent = emptyenv())
    
    vars$zeta <- zeta
    vars$p1 <- p1
    vars$p2 <- p2
    vars$p3 <- p3
    
    if(is.null(cluster)){
      vars$cl <- makeCluster(1)
    } else{
      vars$cl <- cluster
    }
    
    registerDoParallel(vars$cl)                   #making cluster usable with foreach
    vars$no_cores <- length(vars$cl)
    
    #this much will figure out how to interpret the formula
    #first the fixed effects part
    stopifnot(inherits(fixed, "formula"))
    if (missing(data)) {
      barf <- lm(fixed, method = "model.frame")
    } else {
      stopifnot(inherits(data, "data.frame"))
      barf <- lm(fixed, data = data, method = "model.frame")
    }
    x <- model.matrix(fixed, data = barf)
    y <- model.response(barf)
    
    
    #family stuff
    vars$family.glmm<-getFamily(family.glmm)
    
    #check that only binomial has y being matrix
    if(length(dim(y))==2){
      if(vars$family.glmm$family.glmm!="binomial.glmm") {
        stop("For the family you've specified, only a vector is appropriate as the response. Binomial is the only family that allows you to specify the response as a matrix.")
      }	
    }
    
    vars$ntrials <- rep(1, length(y)) #used for Poisson and Bern, essentially untouched
    if(vars$family.glmm$family.glmm=="binomial.glmm"){
      #if the response is a vector, then ntrials stays at 1
      #if the response is a matrix
      if(length(dim(y))==2){
        #make sure it has exactly 2 columns
        if(ncol(y)!=2) stop("Your response must have two columns: the first column reports the number of successes and the second column reports the number of failures.")
        
        # make ntrials a vector with each entry the sum of the entries in the corresponding col of y
        vars$ntrials <- apply(y, MARGIN=1, FUN=sum)
        y <- y[,1] 	#then change y to just be the number of successes
      }
      
      
    }
    #do the check for the specified family
    vars$family.glmm$checkData(y)
    
    
    
    #then the part for the random effects. 
    #first, if it's not a list, make it a list 
    randcall<-random
    if (! is.list(random))
      random <- list(random)
    #put this stuff in a loop and loop along the list
    #for i in 1:length(formula2)
    for (irandom in seq(along = random)) 
    {
      r<-random[[irandom]]
      stopifnot(inherits(r, "formula"))
      if (missing(data)) {
        barf2 <- lm(r, method = "model.frame")
      } else {
        stopifnot(inherits(data, "data.frame"))
        barf2 <- lm(r, data = data, method = "model.frame")
      }
      random[[irandom]] <- model.matrix(r, data = barf2)
      #thisgroup<-varcomps.equal[irandom]
      #names(random)[irandom]<-varcomps.names[thisgroup]
      
      if(length(y)!=nrow(random[[irandom]])) {
        stop("Fixed and random effect model matrices should have same number of rows. This problem sometimes arises due to NAs (missing data).")
      }
    }
    #so now random is a list containing a model matrix for each formula, and some matrices share variance components
    
    #	#family stuff
    #	family.glmm<-getFamily(family.glmm)
    
    #	#check that only binomial has y being matrix
    #	if(length(dim(y))==2){
    #		if(family.glmm$family.glmm!="binomial.glmm") {
    #			stop("For the family you've specified, only a vector is appropriate as the response. Binomial is the only family that allows you to specify the response as a matrix.")
    #		}	
    #	}
    
    #	ntrials <- rep(1, length(y)) #used for Poisson and Bern, essentially untouched
    #	if(family.glmm$family.glmm=="binomial.glmm"){
    #		#if the response is a vector, then ntrials stays at 1
    #		#if the response is a matrix
    #		if(length(dim(y))==2){
    #			#make sure it has exactly 2 columns
    #			if(ncol(y)!=2) stop("Your response must have two columns: the first column reports the number of successes and the second column reports the number of failures.")
    
    #			# make ntrials a vector with each entry the sum of the entries in the corresponding col of y
    #			ntrials <- apply(y, MARGIN=1, FUN=sum)
    #			y <- y[,1] 	#then change y to just be the number of successes
    #		}
    
    
    #	}
    #	#do the check for the specified family
    #	family.glmm$checkData(y)
    
    
    
    if(is.numeric(varcomps.equal)==F) stop("varcomps.equal must be a vector containing numbers to indicate which variance components are equal.")
    if(length(varcomps.equal)!=length(random)){
      stop("The length of varcomps.equal must be equal to the length of the random-effects call.")} 
    if(length(unique(varcomps.equal))!=length(varcomps.names)){
      stop("You must name each unique variance component. Check varcomps.names and varcomps.equal.")} 
    if(min(varcomps.equal)!=1)stop("The vector varcomps.equal must contain numbers starting at 1 to denote which variance components are equal.")	
    levs<-ordered(unique(varcomps.equal))
    
    
    
    
    
    #check p1 p2 p3
    if(!is.numeric(vars$p1))stop("p1 must be a number between 0 and 1")
    if(vars$p1>1) stop("p1 must be a number between 0 and 1")
    if(vars$p1<0) stop("p1 must be a number between 0 and 1")
    if(vars$p1==0) stop("p1 must be nonzero")
    if(!is.numeric(vars$p2))stop("p2 must be a number between 0 and 1")
    if(vars$p2>1) stop("p2 must be a number between 0 and 1")
    if(vars$p2<0) stop("p2 must be a number between 0 and 1")
    if(!is.numeric(vars$p3))stop("p3 must be a number between 0 and 1")
    if(vars$p3>1) stop("p3 must be a number between 0 and 1")
    if(vars$p3<0) stop("p3 must be a number between 0 and 1")
    if(vars$p1+vars$p2+vars$p3!=1) stop("p1+p2+p3 must equal 1")
    
    #this loop is a 2-4-1. We want to check that they're filling in varcomps.equal correctly. 
    #We also want to group all the design matrices that share a variance components.
    #Now z is a list with the number of design mats = number of distinct variance components
    z<-list()
    for(i in 1:length(levs)){
      if(levs[i]!=i) stop("The numbers in the vector varcomps.equal must be consecutive. You must start at 1 and then each entry must be the next consecutive number or a repeat of a previous number.")
      these<-varcomps.equal==i
      thesemats<-random[these]
      z[[i]]<-do.call(cbind,thesemats)
    }
    names(z)<-varcomps.names
    
    if(is.null(weights)){
      wts <- rep(1, length(y))
    } else{
      wts <- weights
    }
    
    #checking weights
    if(typeof(wts) != "double")stop("weights must be a vector")
    vars$wts <- as.vector(wts)
    if(length(vars$wts) != length(y))stop("weights length must match the length of the response vector")
    if(any(!is.numeric(vars$wts)))stop("weights must be a numeric vector")
    if(any(vars$wts < 0))stop("Negative weights not allowed")
    if(any(is.na(vars$wts)))stop("Missing weights not allowed")
    
    #savedx <- x
    #savedy <- y
    #savedw <- weights
    
    #if(any(weights == 0)){
      #drop <- which(weights == 0)
      #if(length(drop) == length(weights))stop("No informative observations. At least one weight must be non-zero")
      #weights <- weights[-drop]
      #x <- x[-drop, , drop=FALSE]
      #y <- if(NCOL(savedy) == 1) y[-drop] else y[-drop, ]
      #for(l in 1:length(z)){
        #z[[l]] <- z[[l]][-drop,] 
      #}
    #}
    
    #w <- sqrt(weights)
    
    #x <- x*w
    #y <- y*w
    
    #for(l in 1:length(z)){
      #for(i in 1:nrow(z[[l]])){
        #for(j in 1:ncol(z[[l]])){
          #z[[l]][i,j] <- z[[l]][i,j]*w[i]
        #}
      #}
    #}
    
    vars$mod.mcml <- list(x=x, y=y, z=z, ntrials=vars$ntrials)
    
    #so now the 3 items are x (matrix), z (list), y (vector)
    #end figuring out how to interpret the formula
    
    #make sure par.init has the right number of parameters and that the variance 
    #components are positive to start. also skip pql bc using par.init instead
    if(!is.null(par.init)){
      nbeta<-ncol(x)
      nbetaplusT<-nbeta+length(z)
      if(length(par.init)!=nbetaplusT) stop("par.init is not the correct length. It should contain initial values for the fixed effects and variance components.")
      vcs<-par.init[-(1:nbeta)]
      if(any(vcs<=10^-9)) stop("Initial values for the variance components in par.init must be positive and sufficiently large (greater than 10^-9).")
      doPQL<-FALSE
      #if par.init is given, we want to use those not PQL
    }
    
    #cache will hold some pql estimates and the importance sampling weights that wouldn't otherwise be returned
    cache <- new.env(parent = emptyenv())
    
    #if the user wants to do pql, do it and use that as the trust start point
    if(doPQL==TRUE){
      #do PQL
      pql.out<-pql(vars$mod.mcml,vars$family.glmm, vars$wts,cache)
      s.pql<-cache$s.twid	
      sigma.pql<-pql.out$sigma
      vars$nu.pql<-sigma.pql^2
      beta.pql<-cache$beta.twid 
      par.init<-c(beta.pql,vars$nu.pql) 
    }
    
    #if the user does not want to do pql, then the best guess of the rand effs is 0
    #if the user did not provide par.init, then an arbitrary guess is 0 for beta
    # and 1 for nu
    if(doPQL==FALSE){
      nrand<-lapply(vars$mod.mcml$z,ncol)
      nrandom<-unlist(nrand)
      totnrandom<-sum(nrandom)
      s.pql<-rep(0,totnrandom)
      if(!is.null(par.init)){ #then par.init is already specified by user
        beta.pql<-par.init[1:nbeta]
        vars$nu.pql<-par.init[-(1:nbeta)]
        sigma.pql<-sqrt(vars$nu.pql)
      }
      if(is.null(par.init)){
        sigma.pql<-vars$nu.pql<-rep(1,length(vars$mod.mcml$z))
        beta.pql<-rep(0,ncol(vars$mod.mcml$x))
        par.init<-c(beta.pql,vars$nu.pql) 
      }
    }
    
    #calculate A*, D* and u*
    nrand<-lapply(vars$mod.mcml$z,ncol)
    nrandom<-unlist(nrand)
    q<-sum(nrandom)
    if(q!=length(s.pql)) stop("Can't happen. Number of random effects returned by PQL must match number of random effects specified by model.")
    eek<-getEk(vars$mod.mcml$z)
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
    
    #now D.star.inv and D.star are both diagonal matrices
    #Diagonal from Matrix package is used bc these are sparse matrices
    #If q (# rand effs) is large, then need to be careful with these
    
    #determine m1, m2, m3 based on probs p1, p2, p3
    foo<-runif(m)
    vars$m1<-sum(foo<vars$p1)
    m2<-sum(foo<vars$p1+vars$p2)-vars$m1	
    m3<-m-vars$m1-m2
    
    #	#generate m1 from N(0,I) and will be scaled to N(0,D) later
    #	zeros<-rep(0,length(u.star))
    #	ones<-rep(1,length(u.star))
    #	ident<-Diagonal(length(u.star),ones)
    #	genData<-genRand(zeros,ident,m1)
    #	
    
    simulate <- function(vars, Dstarnotsparse, m2, m3, beta.pql, D.star.inv){
      #generate m1 from t(0,D*)
      newm1 <- ceiling(vars$m1/vars$no_cores)
      if(vars$m1>0) genData<-rmvt(newm1,sigma=Dstarnotsparse,df=vars$zeta,type=c("shifted"))
      if(vars$m1==0) genData<-NULL		
      
      #generate m2 from N(u*,D*)
      newm2 <- ceiling(m2/vars$no_cores)
      if(m2>0) genData2<-genRand(vars$u.star,vars$D.star,newm2)
      if(m2==0) genData2<-NULL
      
      
      #generate m3 from N(u*,(Z'c''(Xbeta*+zu*)Z+D*^{-1})^-1)
      newm3 <- ceiling(m3/vars$no_cores)
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
        genData3<-genRand(vars$u.star,Sigmuh,newm3)
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
    
    vars$nbeta <- length(beta.pql)
    
    parallel::clusterSetRNGStream(vars$cl, 1234)
    
    clusterEvalQ(vars$cl, library(Matrix))
    clusterEvalQ(vars$cl, library(mvtnorm))
    clusterExport(vars$cl, c("vars", "Dstarnotsparse", "m2", "m3", "beta.pql", "D.star.inv", "simulate"), envir = environment())     #installing variables on each core
    clusterEvalQ(vars$cl, umatparams <- simulate(vars=vars, Dstarnotsparse=Dstarnotsparse, m2=m2, m3=m3, beta.pql=beta.pql, D.star.inv=D.star.inv))
    
    umats <- clusterEvalQ(vars$cl, umatparams$umat)
    umat <- Reduce(rbind, umats)
    vars$newm <- nrow(umat)
    
    #use trust to max the objfun (monte carlo likelihood)
    trust.out<-trust(objfun,parinit=par.init,rinit=10, minimize=FALSE, rmax=rmax, iterlim=iterlim, blather=debug,  cache=cache, vars=vars)
    
    beta.trust<-trust.out$argument[1:length(beta.pql)]
    nu.trust<-trust.out$argument[-(1:length(beta.pql))]
    
    #while(trust.out$converged==FALSE){
    
    #}
    
    names(beta.trust)<-colnames(vars$mod.mcml$x)
    names(nu.trust)<-varcomps.names
    
    if(debug==TRUE){
      debug<-list(beta.pql=beta.pql, nu.pql=vars$nu.pql, D.star=vars$D.star, trust.argpath=trust.out$argpath, u.star=vars$u.star, umat=umat,hessianweights=cache$weights,wtsnumer=cache$numer,wtsdenom=cache$denom,m1=vars$m1,m2=m2,m3=m3,trust.argtry=trust.out$argtry, trust.steptype=trust.out$steptype, trust.accept=trust.out$accept, trust.r=trust.out$r, trust.rho=trust.out$rho, trust.valpath=trust.out$valpath, trust.valtry=trust.out$valtry, trust.preddif=trust.out$preddif, trust.stepnorm=trust.out$stepnorm)
    }
    
    if(is.null(cluster)){
      stopCluster(vars$cl)
    }
    
    return(structure(list(beta=beta.trust,nu=nu.trust, likelihood.value=trust.out$value, likelihood.gradient=trust.out$gradient, likelihood.hessian=trust.out$hessian,
                          trust.converged=trust.out$converged,  mod.mcml=vars$mod.mcml,
                          fixedcall=fixed,randcall=randcall, x=x,y=y, z=random, weights=weights,
                          family.glmm=vars$family.glmm, call=call, varcomps.names=varcomps.names, 
                          varcomps.equal=varcomps.equal, umat=umat, pvec=c(vars$p1, vars$p2, vars$p3), beta.pql=beta.pql, nu.pql=vars$nu.pql, u.pql=vars$u.star, zeta=vars$zeta, cluster=vars$cl, cores=vars$no_cores, debug=debug), class="glmm"))
  }
