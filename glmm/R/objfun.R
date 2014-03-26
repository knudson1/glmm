objfun <-
function(par,nbeta,nu.pql,umat,u.star=u.star,mod.mcml,family.mcml){
	#print(par)
	beta<-par[1:nbeta]
	nu<-par[-(1:nbeta)]
	m<-nrow(umat)

	
	if(sum(nu<=0)>0){
		out<-list(value=-Inf,gradient=rep(1,length(par)),hessian=as.matrix(c(rep(1,length(par)^2)),nrow=length(par)))
	return(out)
	}
	
	Z=do.call(cbind,mod.mcml$z)

	eta<-b<-rep(0,m)
	lfu<-lfu.twid<-lfyu<-list(rep(c(0,0,0),m))
	
	#for each simulated random effect vector
	for(k in 1:m){
		Uk<-umat[k,]  #use the simulated vector as our random effect vec
		eta<-mod.mcml$x%*%beta+Z%*%Uk # calculate eta using it
		zeros<-rep(0,length(Uk))
		lfu[[k]]<-distRand(nu,Uk,mod.mcml$z,zeros)            #log f_theta(u_k)
		lfu.twid[[k]]<-distRand(nu.pql,Uk,mod.mcml$z,u.star)   #log f~_theta(u_k)
		lfyu[[k]]<-el(mod.mcml$y,mod.mcml$x,eta,family.mcml) #log f_theta(y|u_k)
		
		b[k]<-lfu[[k]]$value+lfyu[[k]]$value-lfu.twid[[k]]$value
	}

	a<-max(b)
	thing<-exp(b-a)
	value<-a-log(m)+log(sum(thing))
	
	#weights
	v<-thing/sum(thing)
	
	Gpiece<-matrix(data=NA,nrow=nrow(umat),ncol=length(par))
	
	#lfuky<-NA
	for(k in 1:nrow(umat)){
		Gpiece[k,]<-c(lfyu[[k]]$gradient,lfu[[k]]$gradient)*v[k]	
		
				#lfuky[k]<-c(lfyu[[k]]$gradient,lfu[[k]]$gradient)
		#Gpiece[k,]<-lfuky[k]*v[k]	
	}
	G<-apply(Gpiece,2,sum)
	
	#Hessian has three pieces: panda, lobster, GGT
	panda.list<-list()
	for(k in 1:nrow(umat)){
		panda.list[[k]]<-c(lfyu[[k]]$gradient,lfu[[k]]$gradient)%*%t(c(lfyu[[k]]$gradient,lfu[[k]]$gradient))*v[[k]]

	}
	panda<-addMats(panda.list)
	
	lobster.list<-list()
	for(k in 1:nrow(umat)){

		mat1<-lfyu[[k]]$hessian
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
