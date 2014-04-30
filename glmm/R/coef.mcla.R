coef.mcla <-
function(object,...){
	mod<-object
   	stopifnot(inherits(mod, "mcla"))
	coefficients<-mod$beta
	names(coefficients)<-colnames(mod$x)
	coefficients
}

