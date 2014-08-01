el <-
function(Y,X,eta,family.mcml){
	family.mcml<-getFamily(family.mcml)
	neta<-length(eta)

	if(family.mcml$family.glmm=="bernoulli.glmm"){
		foo<-.C("cum3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(1),cumout=double(neta))$cumout
		mu<-.C("cp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(1),cpout=double(neta))$cpout
		cdub<-.C("cpp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(1),cppout=double(neta))$cppout
	}
	if(family.mcml$family.glmm=="poisson.glmm"){
		foo<-.C("cum3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),cumout=double(neta))$cumout
		mu<-.C("cp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),cpout=double(neta))$cpout
		cdub<-.C("cpp3",eta=as.double(eta),neta=as.integer(neta),type=as.integer(2),cppout=double(neta))$cppout
	}


#	for(i in 1:length(eta)){
#		etai<-eta[i]
#		if(family.mcml$family.glmm=="bernoulli.glmm"){
#			foo[i]<-.C("cum2", eta = as.double(etai), type = as.integer(1),  cumout=double(1))$cumout
#			mu[i]<-.C("cp2",eta=as.double(etai),type=as.integer(1),cpout=double(1))$cpout
#			cdub[i]<-.C("cpp2",eta=as.double(etai),type=as.integer(1),cppout=double(1))$cppout
#		}
#		if(family.mcml$family.glmm=="poisson.glmm"){
#			foo[i]<-.C("cum2", eta = as.double(etai), type = as.integer(2),  cumout=double(1))$cumout
#			mu[i]<-.C("cp2",eta=as.double(etai),type=as.integer(2),cpout=double(1))$cpout
#			cdub[i]<-.C("cpp2",eta=as.double(etai),type=as.integer(2),cppout=double(1))$cppout
#		}
#	}
	value<-Y%*%eta-sum(foo)
	gradient<-t(X)%*%(Y-mu)
	
	cdubmat<-diag(cdub)
	hessian<-t(X)%*%(-cdubmat)%*%X
	
	list(value=value,gradient=gradient,hessian=hessian)
}
