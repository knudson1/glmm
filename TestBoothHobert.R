# test on Booth and Hobert data

library(glmm,lib.loc="glmm.Rcheck")
data(BoothHobert)

#mostly just want to get the PQL estimate 
mod.mcml<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,family.glmm=bernoulli.glmm,m=100) 
beta.pql<-mod.mcml$beta.pql
nu.pql<-mod.mcml$nu.pql
u.pql<-mod.mcml$u.pql
umat<-mod.mcml$umat
x<-mod.mcml$x
y<-mod.mcml$y
z<-mod.mcml$z[[1]]

par<-c(beta.pql,nu.pql)


logfyuk<-function(eta,x,y){
	value<-sum(y*eta)-sum(log(1+exp(eta)))
	Pi<-exp(eta)/(1+exp(eta))
	gradient<-sum(y*x)-sum(x*Pi)
	hessian<-sum(x^2*(Pi-Pi^2)	)
	list(value=value,gradient=gradient,hessian=hessian)
}

#now that I have all the generated rand eff and I know
#their distribution, I can calculate the function that 
#trust is going to maximize
maxthis<-function(par,x,y,z,umat,nu.pql,u.pql){
	beta<-par[1]
	sigsq<-par[2]
	m<-nrow(umat)
	dbb<-db<-b<-rep(0,m)

	#now go through row by row of umat 
	#ie go through each vector of gen rand eff
	for(k in 1:m){
		u<-umat[k,]
		eta<-x*beta+as.vector(z%*%u)
		
		piece1<- logfyuk(eta,x,y)$value	
		piece2<- (-10/2)*(log(sigsq)-sum(u^2)/(2*sigsq))
		piece3<- (-10/2)*(log(nu.pql)-sum((u-u.pql)^2)/(2*nu.pql))	
		b[k]<-piece1+piece2-piece3
	
		#calculate derivs of log f(y|u)
		db[k]<-logfyuk(eta,x,y)$gradient
		dbb[k]<-logfyuk(eta,x,y)$hessian

	}	
	a<-max(b)
	top<-exp(b-a)
	value<-a-log(m)+log(sum(top))
	
	wts<-top/sum(top)
	value
	

}

maxthis(c(mod.mcml$beta,mod.mcml$nu),x=x,y=y,z=z,umat=umat,nu.pql=nu.pql,u.pql=u.pql)
mod.mcml$likelihood.value
