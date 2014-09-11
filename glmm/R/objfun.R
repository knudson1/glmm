objfun <-
function(par,nbeta,nu.pql,umat,u.star=u.star,mod.mcml,family.glmm,cache,distrib,gamm,p1,p2,p3,D.star,Sigmuh){
	beta<-par[1:nbeta]
	nu<-par[-(1:nbeta)]
	m<-nrow(umat)

	if (!missing(cache)) stopifnot(is.environment(cache))
	if(missing(cache)) cache<-new.env(parent = emptyenv())

	if(sum(nu<=0)>0){
		out<-list(value=-Inf,gradient=rep(1,length(par)),hessian=as.matrix(c(rep(1,length(par)^2)),nrow=length(par)))
	return(out)
	}
	
	Z=do.call(cbind,mod.mcml$z)
	T<-length(mod.mcml$z)
	nrand<-lapply(mod.mcml$z,ncol)
	nrandom<-unlist(nrand)


	family.glmm<-getFamily(family.glmm)
	if(family.glmm$family.glmm=="bernoulli.glmm"){family_glmm=1}	
	if(family.glmm$family.glmm=="poisson.glmm"){family_glmm=2}	

	D.star.inv<-solve(D.star)
	Sigmuh.inv<-solve(Sigmuh)
	logdet.D.star.inv<-	sum(log(eigen(D.star.inv,symmetric=TRUE)$values))
	logdet.Sigmuh.inv<-sum(log(eigen(Sigmuh.inv,symmetric=TRUE)$values))
 	myq<-nrow(D.star.inv)

	#for the particular value of nu we're interested in, need to prep for distRandGenC
	eek<-getEk(mod.mcml$z)
	preDinvfornu<-Map("*",eek,(1/nu))
	Dinvfornu<-addVecs(preDinvfornu)
	Dinvfornu<-diag(Dinvfornu)
	logdetDinvfornu<-sum(log(eigen(Dinvfornu,symmetric=TRUE)$values))
	
	meow<-rep(0,T+1)
	meow[1]=0
	meow[2]=nrandom[1]
	if(T>2){	
		for(t in 2:T+1){
		meow[t]=meow[t-1]+nrandom[t-1]
		}
	}
	
	pee<-c(p1,p2,p3)
	n<-nrow(mod.mcml$x)

	stuff<-.C("objfunc", as.double(mod.mcml$y),as.double(t(umat)), as.integer(myq), as.integer(m), as.double(mod.mcml$x), as.integer(n), as.integer(nbeta), as.double(beta), as.double(Z), as.double(Dinvfornu), as.double(logdetDinvfornu),as.integer(family_glmm), as.double(D.star.inv), as.double(logdet.D.star.inv), as.double(u.star), as.double(Sigmuh.inv), as.double(logdet.Sigmuh.inv), pee=as.double(pee), nps=as.integer(length(pee)), T=as.integer(T), nrandom=as.integer(nrandom), meow=as.integer(meow),nu=as.double(nu), v=double(m),value=double(1),gradient=double(length(par)),hessian=double((length(par))^2))

	cache$weights<-stuff$v

	list(value=stuff$value,gradient=stuff$gradient,hessian=matrix(stuff$hessian,ncol=length(par),byrow=FALSE))

}


#objfun <-
#function(par,nbeta,nu.pql,umat,u.star=u.star,mod.mcml,family.glmm,cache,distrib,gamm,p1,p2,p3,D.star,Sigmuh){
#	beta<-par[1:nbeta]
#	nu<-par[-(1:nbeta)]
#	m<-nrow(umat)

#	if (!missing(cache)) stopifnot(is.environment(cache))
#	if(missing(cache)) cache<-new.env(parent = emptyenv())

#	if(sum(nu<=0)>0){
#		out<-list(value=-Inf,gradient=rep(1,length(par)),hessian=as.matrix(c(rep(1,length(par)^2)),nrow=length(par)))
#	return(out)
#	}
#	
#	Z=do.call(cbind,mod.mcml$z)
#	T<-length(mod.mcml$z)
#	nrand<-lapply(mod.mcml$z,ncol)
#	nrandom<-unlist(nrand)

#	eta<-b<-rep(0,m)
#	lfuval<-rep(0,m)
#	lfu.twid<-matrix(data=NA,nrow=m,ncol=4)
#	family.glmm<-getFamily(family.glmm)

#	D.star.inv<-solve(D.star)
#	Sigmuh.inv<-solve(Sigmuh)
#	logdet.D.star.inv<-	sum(log(eigen(D.star.inv,symmetric=TRUE)$values))
#	logdet.Sigmuh.inv<-sum(log(eigen(Sigmuh.inv,symmetric=TRUE)$values))
# 	myq<-nrow(D.star.inv)

#	#for the particular value of nu we're interested in, need to prep for distRandGenC
#	eek<-getEk(mod.mcml$z)
#	preDinvfornu<-Map("*",eek,(1/nu))
#	Dinvfornu<-addVecs(preDinvfornu)
#	Dinvfornu<-diag(Dinvfornu)
#	logdetDinvfornu<-sum(log(eigen(Dinvfornu,symmetric=TRUE)$values))
#	

#	meow<-rep(0,T+1)
#	meow[1]=0
#	meow[2]=nrandom[1]
#	if(T>2){	
#		for(t in 2:T+1){
#		meow[t]=meow[t-1]+nrandom[t-1]
#		}
#	}
#	
#	zeros<-rep(0,ncol(umat))

#	lfyu.val<-rep(0,m)
#	lfyu.gradient<-matrix(0,nrow=m,ncol=length(beta))
#	lfyu.hess<-list()
#	
#	#for each simulated random effect vector
#	for(k in 1:m){
#		Uk<-umat[k,]  #use the simulated vector as our random effect vec	
#		eta<-mod.mcml$x%*%beta+Z%*%Uk # calculate eta using Uk

#		#log f_theta(u_k)
#		lfuval[k]<-.C("distRandGenC",as.double(Dinvfornu),as.double(logdetDinvfornu), as.integer(myq), as.double(Uk), as.double(zeros), double(1))[[6]]

#		#log f_theta(y|u_k) value
#		if(family.glmm$family.glmm=="bernoulli.glmm"){	
#			lfyu.val[k]<-.C("elval",as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta),as.integer(1),value=double(1))$value}
#		if(family.glmm$family.glmm=="poisson.glmm"){		
#			lfyu.val[k]<-.C("elval",as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta),as.integer(2),value=double(1))$value}
#		
#		#value of log f~_theta(u_k)
#		lfu.twid[k,1]<-.C("distRandGenC",as.double(D.star.inv),as.double(logdet.D.star.inv), as.integer(myq), as.double(Uk), as.double(zeros), double(1))[[6]]
#		lfu.twid[k,2]<-.C("distRandGenC",as.double(D.star.inv),as.double(logdet.D.star.inv), as.integer(myq), as.double(Uk), as.double(u.star), double(1))[[6]]
#		lfu.twid[k,3]<-.C("distRandGenC",as.double(Sigmuh.inv),as.double(logdet.Sigmuh.inv), as.integer(myq), as.double(Uk), as.double(u.star), double(1))[[6]]
#		
#		tempmax<-max(lfu.twid[k,1:3])
#		blah<-exp(lfu.twid[k,1:3]-tempmax)
#		pee<-c(p1,p2,p3)
#		qux<-pee%*%blah
#		lfu.twid[k,4]<-tempmax+log(qux)
#		
#		b[k]<-as.numeric(lfuval[k])+as.numeric(lfyu.val[k])-lfu.twid[k,4]
#	}

#	a<-max(b)
#	thing<-exp(b-a)
#	value<-a-log(m)+log(sum(thing))
#	v<-thing/sum(thing)
#	#bk are log weights
#	cache$weights<-exp(b)
#	
#	G<-rep(0,length(par))
#	#Hessian has three pieces: panda, lobster, GGT
#	lobster<-panda<-matrix(rep(0,(length(par))^2),nrow=length(par))

#	for(k in 1:nrow(umat)){
#		Uk<-umat[k,]  #use the simulated vector as our random effect vec	
#		eta<-mod.mcml$x%*%beta+Z%*%Uk # calculate eta using Uk
#		lfutemp<- .C("distRand3C",as.double(nu), as.double(zeros), as.integer(T), as.integer(nrandom), as.integer(meow), as.double(Uk), double(T), double(T^2))  
#		lfugradient<-lfutemp[[7]]
#		lfuhess<-matrix(lfutemp[[8]],byrow=F,nrow=T)

#		#gradient and hessian log f_theta(y|u_k) 
#		if(family.glmm$family.glmm=="bernoulli.glmm"){	
#blah<-.C("elGH",as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta),as.integer(1),gradient=double(ncol(mod.mcml$x)),hessian=double((ncol(mod.mcml$x)^2)))}	
#		if(family.glmm$family.glmm=="poisson.glmm"){		
#blah<-.C("elGH",as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta),as.integer(2),gradient=double(ncol(mod.mcml$x)),hessian=double((ncol(mod.mcml$x)^2)))}
#		lfyu.gradient<-blah$gradient
#		lfyu.hess<-matrix(blah$hessian,byrow=FALSE,nrow=nbeta)

#		#for gradient calculation
#		Gpiece<-c(lfyu.gradient,lfugradient)*v[k]	
#		G<-G+Gpiece

#		#for hessian calculation
#		panda.temp<-c(lfyu.gradient,lfugradient)%*%t(c(lfyu.gradient,lfugradient))*v[k]
#		panda<-panda+panda.temp

#		#for lobster part of hessian calculation
#		newmat<-matrix(data=0,nrow=nbeta+T,ncol=nbeta+T)
#		newmat[1:nbeta,1:nbeta]<-lfyu.hess
#		here<-nbeta+1
#		there<-nbeta+T
#		newmat[here:there,here:there]<-lfuhess	
#		lobster.temp<-newmat*v[k]
#		lobster<-lobster+lobster.temp
#	}
#	
#	#this is for Hessian calculation
#	hessian<-lobster+panda-G%*%t(G)

#	list(value=value,gradient=G,hessian=hessian)
#}

