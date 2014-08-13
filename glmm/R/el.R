el <-
function(Y,X,eta,family.mcml){
	family.mcml<-getFamily(family.mcml)
	neta<-length(eta)

	if(family.mcml$family.glmm=="bernoulli.glmm"){
		foo<-.C("cum3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(1),cumout=double(1))$cumout
		mu<-.C("cp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(1),cpout=double(neta))$cpout
		cdub<-.C("cpp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(1),cppout=double(neta))$cppout
	}
	if(family.mcml$family.glmm=="poisson.glmm"){
		foo<-.C("cum3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),cumout=double(1))$cumout
		mu<-.C("cp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),cpout=double(neta))$cpout
		cdub<-.C("cpp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),cppout=double(neta))$cppout
	}

	#foo2<-.C("matvecmult",a=as.double(Y),b=as.double(eta),nrow=as.integer(1),ncol=as.integer(length(Y)),result=double(1))$result
	foo2<-matvecmult(t(Y),eta)
	value<-foo2-foo
	#value<-as.numeric(Y%*%eta-foo)
	
	gradient<-matvecmult(t(X),Y-mu)
	#gradient<-.C("matvecmult",a=as.double(X),b=as.double(Y-mu),nrow=as.integer(nrow(X)), ncol=as.integer(ncol(X)),result=double(ncol(X)))$result
	#gradient<-t(X)%*%(Y-mu)
	
	cdubmat<-diag(cdub)
	hessian<-t(X)%*%(-cdubmat)%*%X
	
	list(value=value,gradient=gradient,hessian=hessian)
}
