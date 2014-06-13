library(glmm)
data(BoothHobert)
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), data=BoothHobert, family.glmm=bernoulli.glmm, m=100, doPQL=TRUE, debug=TRUE, distrib="normal")

mod.mcml<-mod.mcml1$mod.mcml
z<-mod.mcml$z[[1]]
x<-mod.mcml$x
y<-mod.mcml$y

stuff<-mod.mcml1$debug
beta.pql<-stuff$beta.pql
nu.pql<-stuff$nu.pql
u.star<-stuff$u.star
umat<-stuff$umat
you<-umat[1,]

family.glmm<-bernoulli.glmm

#this should be the same as el
logfyuk<-function(eta,x,y){
	value<-sum(y*eta)-sum(log(1+exp(eta)))
	Pi<-exp(eta)/(1+exp(eta))
	gradient<-sum(y*x)-sum(x*Pi)
	hessian<-sum(x^2*(-Pi+Pi^2)	)
	list(value=value,gradient=gradient,hessian=hessian)
}

#compare el and logfyuk for a value of eta
eta<-rep(2,150)
el<-glmm:::el
this<-el(mod.mcml$y,mod.mcml$x,eta,family.mcml=bernoulli.glmm)	
that<-logfyuk(eta,mod.mcml$x,mod.mcml$y)
all.equal(as.numeric(this$value),as.numeric(that$value))
all.equal(as.numeric(this$gradient),as.numeric(that$gradient))
all.equal(as.numeric(this$hessian),as.numeric(that$hessian))



