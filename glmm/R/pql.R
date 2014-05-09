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
	outer.optim<-suppressWarnings(optim(par=sigma, fn=fn.outer,  beta=beta, s=s, Y=Y ,X=X ,Z=Z ,eek=eek ,family.mcml=family.mcml,cache=cache))

	list(sigma=outer.optim$par)

}


fn.outer <-
function(par,beta,s,Y,X,Z,eek,family.mcml,cache){
	sigma<-par
	Aks<-Map("*",eek,sigma)
	A<-addVecs(Aks) #at this point still a vector
	A<-diag(A) #takes the vector and makes it a diag matrix

	nbeta<-length(beta)
	#run trust
	inner.optim<-trust(fn.inner.trust,parinit=c(beta,s),rinit=5,rmax=10000,minimize=F,Y=Y,X=X,Z=Z,A=A,nbeta=nbeta,
family.mcml=family.mcml,cache=cache)
	
	#get beta and s
	#parms<-inner.optim$argument
	#beta.twid<-cache$beta.twid<-parms[1:nbeta]
	#s.twid<-cache$s.twid<-parms[-(1:nbeta)]
	beta.twid<-cache$beta.twid
	s.twid<-cache$s.twid
	
	#calculate eta.twid
	eta.twid<-(X%*%beta.twid+Z%*%A%*%s.twid)
	
	#calculate el(eta.twid) = piece1
	family.mcml<-getFamily(family.mcml)
	piece1<-el(Y,X,eta.twid,family.mcml)$value
	
	#calculate W = c''(eta.twid)
	dubya<-family.mcml$cpp(eta.twid)
	W<-diag(as.vector(dubya))
	
	#calculate thatthing = A t(Z)WZA+I
	thatthing<-A%*%t(Z)%*%W%*%Z%*%A+diag(nrow(A))
	L<-chol(thatthing)
	
	#calculate piece2a = .5 log det(thatthing)
	piece2a<-sum(log(diag(L)))
	
	
	#calculate piece 2b =  .5 s.twid's.twid 
	piece2b <- .5*s%*%s
	
	#then minvalue= piece2a+piece2b-piece1
	minvalue<-piece2a+piece2b-piece1
	
	as.vector(minvalue)
}

fn.inner.trust <-
function(mypar,Y,X,Z,A,family.mcml,nbeta,cache )
{
	beta<-mypar[1:nbeta]
	s<-mypar[-(1:nbeta)]
	eta<-X%*%beta+Z%*%A%*%s
	family.mcml<-getFamily(family.mcml)
	mu<-family.mcml$cp(eta)

	value<- el(Y,X,eta,family.mcml)$value-.5*s%*%s
	
	#gradient calculation
	db<-t(X)%*%(Y-mu)
	ds<- A%*%t(Z)%*%(Y-mu)-s
	gradient<-c(db,ds)

	#hessian calculation
	cdub<-family.mcml$cpp(eta)
	cdub<-as.vector(cdub)
	cdub<-diag(cdub)
	kyoo<-nrow(A)
	piece1<- (-t(X)%*%cdub%*%X)
	piece2<- (-t(X)%*%cdub%*%Z%*%A)
	piece3<- (-A%*%t(Z)%*%cdub%*%X)
	piece4<- (-A%*%t(Z)%*%cdub%*%Z%*%A -diag(kyoo))
	hessian<- rbind(cbind(piece1,piece2),cbind(piece3,piece4))

	cache$s.twid<-s
	cache$beta.twid<-beta
	
	list(value=value,gradient=gradient,hessian=hessian)
}
