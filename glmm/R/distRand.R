distRand <-
function(nu,U,z.list,mu,distrib){
	# T=number variance components
	T<-length(z.list)
	
	#nrandom is q_t
	nrand<-lapply(z.list,ncol)
	nrandom<-unlist(nrand)
	totnrandom<-sum(nrandom)
	gamma<-totnrandom-1
	
	mu.list<-U.list<-NULL
	if(T==1) {
		U.list[[1]]<-U
		mu.list[[1]]<-mu
		}

	if(T>1){
		U.list[[1]]<-U[1:nrandom[1]] 
		mu.list[[1]]<-mu[1:nrandom[1]]
		for(t in 2:T){
			thing1<-sum(nrandom[1:t-1])+1
			thing2<-sum(nrandom[1:t])
			U.list[[t]]<-U[thing1:thing2]
			mu.list[[t]]<-mu[thing1:thing2]
		}
	}
	
	piece1<-piece2<-piece3<-val<-gradient<-Hessian<-rep(0,T)
	
	#for each variance component
	for(t in 1:T){
		you<-as.vector(U.list[[t]])
		mew<-as.vector(mu.list[[t]])
		Umu<-(you-mew)%*%(you-mew)
		
		if(distrib=="normal"){
			val[t]<- as.numeric(-.5*nrandom[t]*log(nu[t])-Umu/(2*nu[t]))
		
			gradient[t]<- -nrandom[t]/(2*nu[t])+Umu/(2*(nu[t])^2)
		
			Hessian[t]<- nrandom[t]/(2*(nu[t])^2)- Umu/((nu[t])^3)
		}

		if(distrib=="tee"){

			piece1[t]<- as.numeric(-.5*nrandom[t]*log(nu[t]))
			piece2[t]<-as.numeric(Umu/(gamma*nu[t]))
			piece3[t]<- -nrandom[t]/(2*nu[t])
		}

	}
		
	if(distrib=="normal") value<-sum(val)
	if(distrib=="tee") {
		denomstuff<-1+sum(piece2)
		const<-((gamma+totnrandom)/2)
		value<-sum(piece1)+const*log(denomstuff)
		gradient<-piece3+const*piece2/(nu*denomstuff)	
		Hessian<-piece3*nu - const*(piece2*nu/denomstuff)^2-const*(2*piece2/nu^2)/denomstuff
	}
	if(T>1) hessian<-diag(Hessian)
	if(T==1) hessian<-matrix(Hessian,nrow=1,ncol=1)
	
	list(value=value,gradient=gradient,hessian=hessian)		
}
