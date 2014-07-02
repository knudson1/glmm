library(glmm)
data(BoothHobert)
set.seed(1234)
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"), data=BoothHobert, family.glmm=bernoulli.glmm, m=21, doPQL=TRUE, debug=TRUE, distrib="normal")

mod.mcml<-mod.mcml1$mod.mcml
z<-mod.mcml$z[[1]]
x<-mod.mcml$x
y<-mod.mcml$y

stuff<-mod.mcml1$debug
beta.pql<-stuff$beta.pql
nu.pql<-stuff$nu.pql
u.pql<-u.star<-stuff$u.star
umat<-stuff$umat

family.glmm<-bernoulli.glmm

objfun<-glmm:::objfun
getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

############################################
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

############################################
#want to check distRand when we use a normal distribution to get our random effects
distRandCheck<-function(nu,uvec,muvec){
ukmuk<-sum((uvec-muvec)^2)
value<--5*log(nu)-ukmuk/(2*nu)
gradient<- -5/nu +ukmuk/(2*nu^2)
hessian<- 5/(nu^2)-ukmuk/(nu^3)
hessian<-as.matrix(hessian)
list(value=value,gradient=gradient,hessian=hessian)
}

distRand<-glmm:::distRand

you<-umat[1,]
this<-distRandCheck(2,you,u.pql)
that<-distRand(2,you,mod.mcml$z,u.pql)

all.equal(this,that)



###############################################
#will want to check distRandGeneral
#need to get D*
distRandGeneral<-glmm:::distRandGeneral
D.star<-2*diag(10)
this<-distRandGeneral(you,u.pql,D.star)
all.equal(this,that$value)

############################################
#want to check that the value of objfun is the same for a value of nu and beta
m<-nrow(umat)
dbb<-db<-b<-rep(0,m)
sigsq<-nu<-2
beta<-6
Z<-mod.mcml$z[[1]]
D.star.inv<-.5*diag(10)

eta.star<-x*beta.pql+as.vector(Z%*%u.star)
cdouble<-as.vector(bernoulli.glmm()$cpp(eta.star)) #still a vector
cdouble<-diag(cdouble)
Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
Sigmuh<-solve(Sigmuh.inv)

#now go through row by row of umat 
#ie go through each vector of gen rand eff
for(k in 1:m){
	uvec<-umat[k,]
	eta<-x*beta+as.vector(Z%*%uvec)

	piece1<- logfyuk(eta,x,y)$value	
	piece2<- distRandCheck(nu,uvec,rep(0,10))$value
	
	piece3<- distRandGeneral(uvec, rep(0,10),D.star)
	piece4<- distRandGeneral(uvec, u.star, D.star)
	piece5<-distRandGeneral(uvec,u.star,Sigmuh)
	b[k]<-piece1+piece2-1/3*(piece3+piece4+piece5)
	}	
a<-max(b)
top<-exp(b-a)
value<-a+log(mean(top))
#going to compare this against objfun's value
cache<-new.env(parent = emptyenv())
objfun<-glmm:::objfun
that<-objfun(c(beta,nu),nbeta=1,nu.pql=nu.pql,u.star=u.star,mod.mcml=mod.mcml, family.glmm=bernoulli.glmm,cache=cache,distrib="normal",gamm=15,umat=umat,m1=7, m2=7,m3=7,D.star=D.star,Sigmuh=Sigmuh)
all.equal(value,that$value)	
#Given generated random effects, the value of the objective function is correct.
#This plus the test of finite diffs for objfun should be enough.




