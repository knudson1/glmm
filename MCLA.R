require(trust)

##############################################################
#creates list of functions to calculate cumulant for binomial
bernoulli.mcml<-function()
{
	family.mcml<- "bernoulli.mcml"
	cum <- function(eta) log(1+exp(eta))
	cp <- function(eta) exp(eta)/(1+exp(eta))
	cpp<-function(eta) exp(eta)/(1+exp(eta))-exp(2*eta)*(1+exp(eta))^(-2)
	
	checkData<-function(x) {
		bads<-sum(x>1) +sum(x<0)
		if(bads>0) stop("data must be between 0 and 1.")
		return(NULL)	
		}
	out<-list(family.mcml=family.mcml, cum=cum, cp=cp, cpp=cpp, checkData=checkData)
	class(out)<-"mcml.family"
	return(out)
}
##############################################################
#creates list of functions to calculate cumulant for poisson
poisson.mcml<-function()
{
	family.mcml<- "poisson.mcml"
	cum <- function(eta) exp(eta)
	cp <- function(eta) exp(eta)
	cpp<-function(eta) exp(eta)
	checkData<-function(x) {
		bads<-sum(x<0)
		if(bads>0) stop("data must be nonnegative.")
		return(NULL)	
		}
	out<-list(family.mcml=family.mcml, cum=cum,cp=cp,cpp=cpp, checkData=checkData)
	class(out)<-"mcml.family"
	return(out)
}

##############################################################
#figures out what family we're working with. then returns all the 
#stuff associated with that family (such as the cumulant function, 
#its 1st and second derivatives, etc)
getFamily<-function(family.mcml)
{
	if(is.character(family.mcml))
		family.mcml<-get(family.mcml,mode="function",envir=parent.frame())
	if(is.function(family.mcml))
		family.mcml<-family.mcml()
	#right way to check the class if there might be more than one class (maybe bc of hierarchies)
	if(!inherits(family.mcml,"mcml.family")) 
		stop(" 'family' not recognized") 
	
	family.mcml	
}
##############################################################
# calculates the log of the data density. l(eta), its gradient
# vector, and the hessian matrix 
el<-function(Y,X,eta,family.mcml){
	family.mcml<-getFamily(family.mcml)
	value<-Y%*%eta-sum(family.mcml$cum(eta))
	
	mu<-family.mcml$cp(eta)
	gradient<-t(X)%*%(Y-mu)
	
	cdub<-as.vector(family.mcml$cpp(eta))
	cdubmat<-diag(cdub)
	hessian<-t(X)%*%(-cdubmat)%*%X
	
	list(value=value,gradient=gradient,hessian=hessian)
}


##############################################################
# concatenates vectors of same length. works for any number of vectors. each vector is an item in a list. returns the summed vector. sums component-wise
#created this for summing diagonals of diagonal matrices, but works for any vectors of uniform length
#this has been checked by doing simple examples to make sure it works properly.
catVecs<-function(vecs){
	 if (! is.list(vecs))
        vecs <- list(vecs)

	lengths <-lapply(vecs,length)
	lengths<-unlist(lengths)
	out<-rep(0,sum(lengths))
	starthere<-1
	
	for(d in 1:length(vecs)){
		endhere<-starthere+lengths[d]-1
		out[c(starthere:endhere)]<-unlist(vecs[[d]])
		starthere<-starthere+lengths[d]
	}
	out
}


##############################################################
# adds vectors of same length. works for any number of vectors. each vector is an item in a list. returns the summed vector. sums component-wise
#created this for summing diagonals of diagonal matrices, but works for any vectors of uniform length
#this has been checked by doing simple examples to make sure it works properly.
addVecs<-function(vecs){
	 if (! is.list(vecs))
        vecs <- list(vecs)

	dvecs <-length(vecs[[1]])
	out<-rep(0,dvecs)
	
	
	for(d in 1:dvecs){
		thing<-lapply(vecs,"[[",d)
		thing<-unlist(thing)
		out[d]<-sum(thing)
	}
	out
}

##############################################################
addMats<-function(matList){
		 if (! is.list(matList))
        matList <- list(matList)
        
        ncells<-length(as.vector(matList[[1]]))
        nr<-nrow(matList[[1]])
        matAsVec<-lapply(matList,as.vector)
        
        out<-rep(0,ncells)

        for(d in 1:ncells){
        	thing<-lapply(matAsVec,"[[",d)
        	thing<-unlist(thing)
        	out[d]<-sum(thing)
        }
	matrix(data=out,nrow=nr)
}


##############################################################
# the log density of the random effects, its gradient vector, and 
# the hessian matrix.
# Finite differences check works.
distRand<-function(nu,U,z.list,mu){
	# T=number variance components
	T<-length(z.list)
	
	#nrandom is q_t
	nrand<-lapply(z.list,ncol)
	nrandom<-unlist(nrand)
	totnrandom<-sum(nrandom)
	
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
	
	val<-gradient<-Hessian<-rep(0,T)
	
	#for each variance component
	for(t in 1:T){
		you<-as.vector(U.list[[t]])
		mew<-as.vector(mu.list[[t]])
		Umu<-(you-mew)%*%(you-mew)
		val[t]<- as.numeric(-.5*nrandom[t]*log(nu[t])-Umu/(2*nu[t]))
		
		gradient[t]<- -nrandom[t]/(2*nu[t])+Umu/(2*(nu[t])^2)
		
		Hessian[t]<- nrandom[t]/(2*(nu[t])^2)- Umu/((nu[t])^3)
		
	}
		
	value<-sum(val)
	if(T>1) hessian<-diag(Hessian)
	if(T==1) hessian<-matrix(Hessian,nrow=1,ncol=1)
	
	list(value=value,gradient=gradient,hessian=hessian)		
}



##############################################################
# when doing PQL, we have an inner optimization and outer optimization.
# this is the inner optimization's function value. this will be MAXIMIZED for trust.
fn.inner.trust<-function(mypar,Y,X,Z,A,family.mcml,nbeta )
{
	beta<-mypar[1:nbeta]
	s<-mypar[-(1:nbeta)]
	eta<-X%*%beta+Z%*%A%*%s
	family.mcml<-getFamily(family.mcml)
	mu<-family.mcml$cp(eta)

	value<- el(Y,X,eta,family.mcml)$value-.5*s%*%s
	
	#gradient calculation
	db<-t(X)%*%(Y-mu)
	ds<- A%*%t(Z)%*%(Y-mu)-s
	gradient<-c(db,ds)

	#hessian calculation
	cdub<-family.mcml$cpp(eta)
	cdub<-as.vector(cdub)
	cdub<-diag(cdub)
	kyoo<-nrow(A)
	piece1<- (-t(X)%*%cdub%*%X)
	piece2<- (-t(X)%*%cdub%*%Z%*%A)
	piece3<- (-A%*%t(Z)%*%cdub%*%X)
	piece4<- (-A%*%t(Z)%*%cdub%*%Z%*%A -diag(kyoo))
	hessian<- rbind(cbind(piece1,piece2),cbind(piece3,piece4))
	
	list(value=value,gradient=gradient,hessian=hessian)
}






##############################################################
# when doing PQL, we have an inner optimization and outer optimization.
# this is the outer optimization's function value

fn.outer<-function(par,beta,s,Y,X,Z,eek,family.mcml){
	sigma<-par
	Aks<-Map("*",eek,sigma)
	A<-addVecs(Aks) #at this point still a vector
	A<-diag(A) #takes the vector and makes it a diag matrix

	nbeta<-length(beta)
	#run trust
	inner.optim<-trust(fn.inner.trust,parinit=c(beta,s),rinit=5,rmax=10000,minimize=F,Y=Y,X=X,Z=Z,A=A,nbeta=nbeta,family.mcml=family.mcml)
	
	#get beta and s
	parms<-inner.optim$argument
	beta.twid<<-parms[1:nbeta]
	s.twid<<-parms[-(1:nbeta)]
	
	#calculate eta.twid
	eta.twid<-(X%*%beta.twid+Z%*%A%*%s.twid)
	
	#calculate el(eta.twid) = piece1
	family.mcml<-getFamily(family.mcml)
	piece1<-el(Y,X,eta.twid,family.mcml)$value
	
	#calculate W = c''(eta.twid)
	dubya<-family.mcml$cpp(eta.twid)
	W<-diag(as.vector(dubya))
	
	#calculate thatthing = A t(Z)WZA+I
	thatthing<-A%*%t(Z)%*%W%*%Z%*%A+diag(nrow(A))
	L<-chol(thatthing)
	
	#calculate piece2a = .5 log det(thatthing)
	piece2a<-sum(log(diag(L)))
	
	
	#calculate piece 2b =  .5 s.twid's.twid 
	piece2b <- .5*s%*%s
	
	#then minvalue= piece2a+piece2b-piece1
	minvalue<-piece2a+piece2b-piece1
	
	as.vector(minvalue)
}

##############################################################
# PQL
pql<-function(mod.mcml,family.mcml){
	eek<-getEk(mod.mcml$z)
	
	#need inits for parameters
	sigma<-rep(1,length(mod.mcml$z))
	beta<-rep(0,ncol(mod.mcml$x))
	
	#need inits for random effects
	nrand<-lapply(mod.mcml$z,ncol)
	nrandom<-unlist(nrand)
	totnrandom<-sum(nrandom)
	s<-rep(1,totnrandom)
	
	#need to give this Y, X, etc
	Y=mod.mcml$y
	X=mod.mcml$x
	Z=do.call(cbind,mod.mcml$z)
	
	#outer.optim<-optim(par=sigma, fn=fn.outer,  beta=beta, s=s, Y=Y ,X=X ,Z=Z ,eek=eek ,family.mcml=family.mcml)
	outer.optim<-suppressWarnings(optim(par=sigma, fn=fn.outer,  beta=beta, s=s, Y=Y ,X=X ,Z=Z ,eek=eek ,family.mcml=family.mcml))

	list(sigma=outer.optim$par,beta=beta.twid,s=s.twid)

}


##############################################################
#gets the diagonal of the Ek matrices. each Ek diagonal is an item in the list.
# Recall that A= sum_k Ek*sigma, where A^2 = D=var(u)
#The number of items in the list Ek = the number of variance components
#The length of each diagonal vector is equal to the number of random effects.
#input: list z containing model matrices 
#each z has n rows. ncols of z is the number of random effects
#we see this in mod.mcml$z
#output: the list of Ek diagonals.

getEk<-function(z){
	nrand<-lapply(z,ncol)
	nrandom<-unlist(nrand)
	totnrandom<-sum(nrandom)
	nsigma<-length(z)
	
	Ekmat<-matrix(data=c(0),nrow=nsigma,ncol=totnrandom)
	Ek<-list()

	cumnrand<-c()
	cumnrand[1]<-0
	for(i in 1:length(nrandom)){
		cumnrand[i+1]<-sum(nrandom[1:i])
	}
	
	#each row of Ekmat is a diagonal of one of the Eks

	for(i in 1:nsigma){
		for(j in 1:totnrandom){
			if(cumnrand[i]<j&cumnrand[i+1]>j-1)		Ekmat[i,j]<-1
		}
		Ek[[i]]<-Ekmat[i,]
	}
	Ek
}




##############################################################
#genRand generates the random effects based on the estimate of the variance component (nuEst)
#sigma.star is from pql
#s.star is from pql (standardized random effects from pql)
#z.list is the list of model matrices (one for each var camp)
#m is the monte carlo sample size
#returns u (generated random effects) and u.star (mean of the u's)
genRand<-function(sigma.star,s.star,z.list,m){
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




##############################################################
#to be maximized by trust
#finite differences show all three pieces are ok!
objfun<-function(par,nbeta,nu.pql,umat,u.star=u.star,mod.mcml,family.mcml){
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






##############################################################
#varcomps.names are the names of the variance components
#varcomps.equal is an optional argument
mcml<-function(fixed,random,varcomps.names,data,family.mcml,m,varcomps.equal,doPQL=T){
	if(missing(varcomps.names)) stop("Names for the variance components must be supplied through varcomps.names")
	if(is.vector(varcomps.names)!=1) stop("varcomps.names must be a vector")

	if(missing(varcomps.equal)){
		varcomps.equal<- c(1:length(varcomps.names))}
	call<-match.call()

	#this much will figure out how to interpret the formula
	#first the fixed effects part
	stopifnot(inherits(fixed, "formula"))
    if (missing(data)) {
        barf <- lm(fixed, method = "model.frame")
    } else {
        stopifnot(inherits(data, "data.frame"))
        barf <- lm(fixed, data = data, method = "model.frame")
    }
    x <- model.matrix(fixed, data = barf)
    y <- model.response(barf)
    
    #then the part for the random effects. 
    #first, if it's not a list, make it a list 
    randcall<-random
	if (! is.list(random))
        random <- list(random)
    #put this stuff in a loop and loop along the list
    #for i in 1:length(formula2)
    for (irandom in seq(along = random)) {
    		r<-random[[irandom]]
    		stopifnot(inherits(r, "formula"))
    		if (missing(data)) {
       		 barf2 <- lm(r, method = "model.frame")
   		} else {
       		 stopifnot(inherits(data, "data.frame"))
        	barf2 <- lm(r, data = data, method = "model.frame")
    		}
    	random[[irandom]] <- model.matrix(r, data = barf2)
	#thisgroup<-varcomps.equal[irandom]
	#names(random)[irandom]<-varcomps.names[thisgroup]

	if(length(y)!=nrow(random[[irandom]])) stop("Fixed and random effect model matrices should have same number of rows")
	}
	#so now random is a list containing a model matrix for each formula, and some matrices share variance components

	if(is.numeric(varcomps.equal)==F) stop("varcomps.equal must be a vector containing numbers to indicate which variance components are equal.")
	if(length(varcomps.equal)!=length(random)){
		stop("The length of varcomps.equal must be equal to the length of the random-effects call.")} 
	if(length(unique(varcomps.equal))!=length(varcomps.names)){
		stop("You must name each unique variance component. Check varcomps.names and varcomps.equal.")} 
	if(min(varcomps.equal)!=1)stop("The vector varcomps.equal must contain numbers starting at 1 to denote which variance components are equal.")	
	levs<-ordered(unique(varcomps.equal))
	
	#this loop is a 2-4-1. We want to check that they're filling in varcomps.equal correctly. 
	#We also want to group all the design matrices that share a variance components.
	#Now z is a list with the number of design mats = number of distinct variance components
	z<-list()
	for(i in 1:length(levs)){
		if(levs[i]!=i) stop("The numbers in the vector varcomps.equal must be consecutive. You must start at 1 and then each entry must be the next consecutive number or a repeat of a previous number.")
		these<-varcomps.equal==i
		thesemats<-random[these]
		z[[i]]<-do.call(cbind,thesemats)
	}
	names(z)<-varcomps.names

	mod.mcml<-structure(list(x = x, z=z,y = y), class = "bar")
   	mod.mcml
	
	#so now the 3 items are x (matrix), z (list), y (vector)
	#end figuring out how to interpret the formula
	
	if(doPQL==T){
	      #do PQL
	      pql.out<-pql(mod.mcml,family.mcml)
	      s.pql<-pql.out$s	
	      sigma.pql<-pql.out$sigma
	      nu.pql<-sigma.pql^2
	      beta.pql<-pql.out$beta

	}
	
	if(doPQL==F){
	      nrand<-lapply(mod.mcml$z,ncol)
	      nrandom<-unlist(nrand)
	      totnrandom<-sum(nrandom)
	      s.pql<-rep(0,totnrandom)
	      nu.pql<-rep(1,length(mod.mcml$z))
	      beta.pql<-rep(1,ncol(mod.mcml$x))
	}
	
	par.init<-c(pql.out$beta,nu.pql) 
	
	# generate random effects
	genData<-genRand(sigma.pql,s.pql,mod.mcml$z,m)
	umat<-genData$u
	u.star<-genData$u.star
	
	#use trust to max the objfun (monte carlo likelihood)
	trust.out<-trust(objfun,parinit=par.init,rinit=10, rmax=10000, 
iterlim=100, minimize=F, nbeta=length(beta.pql), nu.pql=nu.pql, 
umat=umat, mod.mcml=mod.mcml, family.mcml=family.mcml, m=m,u.star=u.star,blather=T)
	
	beta.trust<-trust.out$argument[1:length(beta.pql)]
	nu.trust<-trust.out$argument[-(1:length(beta.pql))]

	trust.argpath<-trust.out$argpath

	names(beta.trust)<-colnames(mod.mcml$x)
	names(nu.trust)<-varcomps.names
	
	return(structure(list(beta=beta.trust,nu=nu.trust, likelihood.value=trust.out$value, 
likelihood.gradient=trust.out$gradient, likelihood.hessian=trust.out$hessian,trust.converged=trust.out$converged, beta.pql=beta.pql, nu.pql=nu.pql, mod.mcml=mod.mcml, trust.argpath=trust.argpath, fixedcall=fixed,randcall=randcall, x=x,y=y, z=random, family.mcml=family.mcml, call=call,umat=umat, varcomps.names=varcomps.names, varcomps.equal=varcomps.equal,u.pql=u.star), class="mcla"))
}



##############################################################
summary.mcla<-function(mod.mcml){
    stopifnot(inherits(mod.mcml, "mcla"))

	fixedcall<-mod.mcml$fixedcall
	randcall<-mod.mcml$randcall
	call<-mod.mcml$call
	x<-mod.mcml$x
	y<-mod.mcml$y
	z<-mod.mcml$z
	
	#the coefficients matrix for fixed effects
	beta<-mod.mcml$beta
	nbeta<-length(beta)
	hessian<-mod.mcml$likelihood.hessian
	if(det(hessian)==0) {
            warning(paste("estimated Fisher information matrix not positive",
               "definite, making all standard errors infinite"))
            all.ses <- rep(Inf, nrow(hessian))
        }
	else{	all.ses<-sqrt(diag(solve(-hessian)))}
	
	beta.se<-all.ses[1:nbeta]
	zval<-beta/beta.se
	coefmat<-cbind(beta,beta.se,zval,2*pnorm(abs(zval),lower.tail=F))
	colnames(coefmat)<-c("Estimate","Std. Error", "z value", "Pr(>|z|)")
	rownames(coefmat)<-colnames(mod.mcml$x)
	
	nu<-mod.mcml$nu
	nu.se<-all.ses[-(1:nbeta)]
	nuzval<-nu/nu.se
	nucoefmat<-cbind(nu,nu.se,nuzval,pnorm(abs(nuzval),lower.tail=F))
	colnames(nucoefmat)<-c("Estimate","Std. Error", "z value", "Pr(>|z|)/2")
	rownames(nucoefmat)<-mod.mcml$varcomps.names
	
	return(structure(list(x=x,y=y, z=z, coefmat=coefmat, fixedcall=fixedcall, randcall=randcall,
family.mcml=mod.mcml$family.mcml,call=call,nucoefmat=nucoefmat),class="summary.mcla"))
}

##############################################################
print.summary.mcla<-function(summ){
    stopifnot(inherits(summ, "summary.mcla"))

    cat("\nCall:\n", paste(deparse(summ$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
   
cat("Fixed Effects:")
   cat("\n")

	printCoefmat(summ$coefmat,digits=3)
   cat("\n")

   cat("\n")
cat("Variance Components for Random Effects (P-values are one-tailed):")
   cat("\n")

	printCoefmat(summ$nucoefmat,digits=3)
   cat("\n")

}






