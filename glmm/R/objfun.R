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

	eta<-b<-rep(0,m)
	lfu<-lfyu<-list(rep(c(0,0,0),m))
	lfu.twid<-matrix(data=NA,nrow=m,ncol=4)
	family.glmm<-getFamily(family.glmm)
	
	#for each simulated random effect vector
	for(k in 1:m){
		Uk<-umat[k,]  #use the simulated vector as our random effect vec	
		eta<-mod.mcml$x%*%beta+Z%*%Uk # calculate eta using it
		zeros<-rep(0,length(Uk))

		#log f_theta(u_k)
		lfu[[k]]<-distRand(nu,Uk,mod.mcml$z,zeros) 

		#log f_theta(y|u_k)
		if(family.glmm$family.glmm=="bernoulli.glmm"){	
			blah<-.C("elc",as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta),as.integer(1),value=double(1),gradient=double(ncol(mod.mcml$x)),hessian=double((ncol(mod.mcml$x)^2)))}
		if(family.glmm$family.glmm=="poisson.glmm"){		
			blah<-.C("elc",as.double(mod.mcml$y),as.double(mod.mcml$x),as.integer(nrow(mod.mcml$x)),as.integer(ncol(mod.mcml$x)),as.double(eta),as.integer(2),value=double(1),gradient=double(ncol(mod.mcml$x)),hessian=double((ncol(mod.mcml$x)^2)))}
		lfyu[[k]]<-list(blah$value,blah$gradient,matrix(blah$hessian,byrow=FALSE,nrow=nbeta))
		#log f~_theta(u_k)
		if(distrib=="normal"){
			lfu.twid[k,1]<-distRandGeneral(Uk,zeros,D.star)}
		if(distrib=="tee") {
			lfu.twid[k,1]<-tdist(nu.pql,Uk,mod.mcml$z,u.star,gamm)}
		lfu.twid[k,2]<-distRandGeneral(Uk,u.star,D.star)
		lfu.twid[k,3]<-distRandGeneral(Uk,u.star,Sigmuh)
		
		tempmax<-max(lfu.twid[k,1:3])
		blah<-exp(lfu.twid[k,1:3]-tempmax)
		pee<-c(p1,p2,p3)
		qux<-pee%*%blah
		lfu.twid[k,4]<-tempmax+log(qux)
		
		b[k]<-as.numeric(lfu[[k]]$value)+as.numeric(lfyu[[k]][[1]])-lfu.twid[k,4]
	}

	a<-max(b)
	thing<-exp(b-a)
	value<-a-log(m)+log(sum(thing))
	v<-thing/sum(thing)
	#bk are log weights
	cache$weights<-exp(b)
	
	Gpiece<-matrix(data=NA,nrow=nrow(umat),ncol=length(par))
	
	#lfuky<-NA
	for(k in 1:nrow(umat)){
		Gpiece[k,]<-c(lfyu[[k]][[2]],lfu[[k]]$gradient)*v[k]	
		
				#lfuky[k]<-c(lfyu[[k]]$gradient,lfu[[k]]$gradient)
		#Gpiece[k,]<-lfuky[k]*v[k]	
	}
	G<-apply(Gpiece,2,sum)
	
	#Hessian has three pieces: panda, lobster, GGT
	panda.list<-list()
	for(k in 1:nrow(umat)){
		panda.list[[k]]<-c(lfyu[[k]][[2]],lfu[[k]]$gradient)%*%t(c(lfyu[[k]][[2]],lfu[[k]]$gradient))*v[[k]]

	}
	panda<-addMats(panda.list)
	
	lobster.list<-list()
	for(k in 1:nrow(umat)){

		mat1<-lfyu[[k]][[3]]
		mat2<-lfu[[k]]$hessian

		d1<-nrow(mat1)
		d2<-nrow(mat2)
		newmat<-matrix(data=0,nrow=d1+d2,ncol=d1+d2)

		newmat[1:d1,1:d1]<-mat1
		here<-d1+1
		there<-d1+d2
		newmat[here:there,here:there]<-mat2	
		lobster.list[[k]]<-newmat*v[k]
	}
	lobster<-addMats(lobster.list)
	
	hessian<-lobster+panda-G%*%t(G)
	list(value=value,gradient=G,hessian=hessian)
}
