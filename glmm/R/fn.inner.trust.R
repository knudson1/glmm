fn.inner.trust <-
function(mypar,Y,X,Z,A,family.mcml,nbeta )
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
	
	list(value=value,gradient=gradient,hessian=hessian)
}
