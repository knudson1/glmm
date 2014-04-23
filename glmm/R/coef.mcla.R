coef.mcla <-
function(mod){
   	stopifnot(inherits(mod, "mcla"))
	coefficients<-mod$beta
	names(coefficients)<-colnames(mod$x)
	coefficients
}
