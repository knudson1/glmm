summary.mcla <-
function(object){
    mod.mcml<-object
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
	coefficients<-coefmat[,1]
	
	nu<-mod.mcml$nu
	nu.se<-all.ses[-(1:nbeta)]
	nuzval<-nu/nu.se
	nucoefmat<-cbind(nu,nu.se,nuzval,pnorm(abs(nuzval),lower.tail=F))
	colnames(nucoefmat)<-c("Estimate","Std. Error", "z value", "Pr(>|z|)/2")
	rownames(nucoefmat)<-mod.mcml$varcomps.names
	
	return(structure(list(x=x,y=y, z=z, coefmat=coefmat, fixedcall=fixedcall, randcall=randcall, coefficients=coefficients,
family.mcml=mod.mcml$family.mcml,call=call,nucoefmat=nucoefmat),class="summary.mcla"))
}
