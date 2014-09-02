#all normal distRand versions now moved to C and to the test files
tdist <-function(nu,U,z.list,mu,gamm){
	#use nu and z.list to get D which is scale matrix
	eek<-getEk(z.list)
	Dvecs<-Map("*",eek,nu)
	Dvec<-addVecs(Dvecs) #at this point still a vector
	Dmat<-diag(Dvec)

	value<-dmvt(U,delta=mu,sigma=Dmat,log=TRUE,df=gamm,type="shifted")
	list(value=value)

}

# no longer used
#distRandGeneral2<-function(uvec,mu,Sigma.inv,logDetSigmaInv){
#	print("still using distRandGeneral2.R")
#	umu<-uvec-mu
#	piece2<-t(umu)%*%Sigma.inv%*%umu

#	as.vector(.5*(logDetSigmaInv-piece2))
#}

# old version of distRand, used only for tests




