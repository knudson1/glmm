genRand <-
function(sigma.star,s.star,z.list,m){
	nrand<-lapply(z.list,ncol)
	nrandom<-unlist(nrand)
	q<-sum(nrandom)
	if(q!=length(s.star)) stop("Number of random effects in dispute between s.star and mod.mcml$z")
	
	eek<-getEk(z.list)
	Aks<-Map("*",eek,sigma.star)
	A.star<-addVecs(Aks) #at this point still a vector
	#A.star<-diag(A.star)
	u.star<-A.star*s.star
	
	u<-matrix(data=NA,nrow=m,ncol=q)
	for(k in 1:m){
		u.swoop<-rnorm(q)
		u[k,]<-u.swoop*A.star+u.star
	}
	
	list(u=u,u.star=u.star)
}
