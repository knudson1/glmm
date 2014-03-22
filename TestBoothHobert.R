# test on Booth and Hobert data

source("MCLA.R")
dat<-read.csv("BoothHobertData.csv")
 dat$z1<-as.factor(dat$z1)

#mostly just want to get the PQL estimate 
mod.mcml<-mcml(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=dat,family.mcml=bernoulli.mcml,m=4) 
beta.pql<-mod.mcml$beta.pql
nu.pql<-mod.mcml$nu.pql
u.pql<-mod.mcml$u.pql
umat<-mod.mcml$umat
x<-mod.mcml$x
y<-mod.mcml$y
z<-mod.mcml$z[[1]]

par<-c(beta.pql,nu.pql)

#now that I have all the generated rand eff and I know
#their distribution, I can calculate the function that 
#trust is going to maximize
maxthis<-function(par,x,y,z,umat,nu.pql,u.pql){
	beta<-par[1]
	sigsq<-par[2]
	m<-nrow(umat)
	b<-rep(0,m)

	#now go through row by row of umat 
	#ie go through each vector of gen rand eff
	for(i in 1:m){
		u<-umat[i,]
		eta<-x*beta+as.vector(z%*%u)
		
		piece1<- y%*%exp(eta)-sum(log(1+exp(eta)))	
		piece2<- (-10/2)*(log(sigsq)-sum(u^2)/(2*sigsq))
		piece3<- (-10/2)*(log(nu.pql)-sum((u-u.pql)^2)/(2*nu.pql))	
		b[i]<-piece1+piece2-piece3
	}	
	a<-max(b)
	value<-a-log(m)+log(sum(exp(b-a)))
	value
}