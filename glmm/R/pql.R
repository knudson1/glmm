pql <-
function(mod.mcml,family.mcml,cache){
   	if (! missing(cache))
       	    stopifnot(is.environment(cache))
	eek<-getEk(mod.mcml$z)
	
	#need inits for parameters
	sigma<-rep(1,length(mod.mcml$z))
	beta<-rep(0,ncol(mod.mcml$x))
	
	#need inits for random effects
	nrand<-lapply(mod.mcml$z,ncol)
	nrandom<-unlist(nrand)
	totnrandom<-sum(nrandom)
	s<-rep(1,totnrandom)
	
	#need to give this Y, X, etc
	Y=mod.mcml$y
	X=mod.mcml$x
	Z=do.call(cbind,mod.mcml$z)
	
	#outer.optim<-optim(par=sigma, fn=fn.outer,  beta=beta, s=s, Y=Y ,X=X ,Z=Z ,eek=eek ,family.mcml=family.mcml,cache=cache)
	outer.optim<-suppressWarnings(optim(par=sigma, fn=fn.outer,  beta=beta, s=s, Y=Y ,X=X ,Z=Z ,eek=eek ,family.mcml=family.mcml))

	list(sigma=outer.optim$par)

}
