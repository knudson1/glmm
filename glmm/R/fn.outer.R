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
