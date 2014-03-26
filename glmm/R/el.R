el <-
function(Y,X,eta,family.mcml){
	family.mcml<-getFamily(family.mcml)
	value<-Y%*%eta-sum(family.mcml$cum(eta))
	
	mu<-family.mcml$cp(eta)
	gradient<-t(X)%*%(Y-mu)
	
	cdub<-as.vector(family.mcml$cpp(eta))
	cdubmat<-diag(cdub)
	hessian<-t(X)%*%(-cdubmat)%*%X
	
	list(value=value,gradient=gradient,hessian=hessian)
}
