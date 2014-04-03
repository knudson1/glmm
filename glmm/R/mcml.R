mcml <-
function(fixed,random,varcomps.names,data,family.mcml,m,varcomps.equal,doPQL=T){
	if(missing(varcomps.names)) stop("Names for the variance components must be supplied through varcomps.names")
	if(is.vector(varcomps.names)!=1) stop("varcomps.names must be a vector")

	if(missing(varcomps.equal)){
		varcomps.equal<- c(1:length(varcomps.names))}
	call<-match.call()

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
    
    #then the part for the random effects. 
    #first, if it's not a list, make it a list 
    randcall<-random
	if (! is.list(random))
        random <- list(random)
    #put this stuff in a loop and loop along the list
    #for i in 1:length(formula2)
    for (irandom in seq(along = random)) {
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

	if(length(y)!=nrow(random[[irandom]])) stop("Fixed and random effect model matrices should have same number of rows")
	}
	#so now random is a list containing a model matrix for each formula, and some matrices share variance components

	if(is.numeric(varcomps.equal)==F) stop("varcomps.equal must be a vector containing numbers to indicate which variance components are equal.")
	if(length(varcomps.equal)!=length(random)){
		stop("The length of varcomps.equal must be equal to the length of the random-effects call.")} 
	if(length(unique(varcomps.equal))!=length(varcomps.names)){
		stop("You must name each unique variance component. Check varcomps.names and varcomps.equal.")} 
	if(min(varcomps.equal)!=1)stop("The vector varcomps.equal must contain numbers starting at 1 to denote which variance components are equal.")	
	levs<-ordered(unique(varcomps.equal))
	
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

	mod.mcml<-structure(list(x = x, z=z,y = y), class = "bar")
   	mod.mcml
	
	#so now the 3 items are x (matrix), z (list), y (vector)
	#end figuring out how to interpret the formula
	
	if(doPQL==T){
	      #do PQL
	      pql.out<-pql(mod.mcml,family.mcml)
	      s.pql<-pql.out$s	
	      sigma.pql<-pql.out$sigma
	      nu.pql<-sigma.pql^2
	      beta.pql<-pql.out$beta

	}
	
	if(doPQL==F){
	      nrand<-lapply(mod.mcml$z,ncol)
	      nrandom<-unlist(nrand)
	      totnrandom<-sum(nrandom)
	      s.pql<-rep(0,totnrandom)
	      nu.pql<-rep(1,length(mod.mcml$z))
	      beta.pql<-rep(1,ncol(mod.mcml$x))
	}
	
	par.init<-c(pql.out$beta,nu.pql) 
	
	# generate random effects
	genData<-genRand(sigma.pql,s.pql,mod.mcml$z,m)
	umat<-genData$u
	u.star<-genData$u.star
	
	#use trust to max the objfun (monte carlo likelihood)
	trust.out<-trust(objfun,parinit=par.init,rinit=10, rmax=10000, 
iterlim=100, minimize=F, nbeta=length(beta.pql), nu.pql=nu.pql, 
umat=umat, mod.mcml=mod.mcml, family.mcml=family.mcml, m=m,u.star=u.star,blather=T)
	
	beta.trust<-trust.out$argument[1:length(beta.pql)]
	nu.trust<-trust.out$argument[-(1:length(beta.pql))]

	trust.argpath<-trust.out$argpath

	names(beta.trust)<-colnames(mod.mcml$x)
	names(nu.trust)<-varcomps.names
	
	return(structure(list(beta=beta.trust,nu=nu.trust, likelihood.value=trust.out$value, 
likelihood.gradient=trust.out$gradient, likelihood.hessian=trust.out$hessian,
	trust.converged=trust.out$converged, beta.pql=beta.pql, nu.pql=nu.pql, mod.mcml=mod.mcml,
	trust.argpath=trust.argpath, fixedcall=fixed,randcall=randcall, x=x,y=y, z=random,
	family.mcml=family.mcml, call=call,umat=umat, varcomps.names=varcomps.names, 
	varcomps.equal=varcomps.equal,u.pql=u.star), class="mcla"))
}
