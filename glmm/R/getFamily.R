getFamily <-function(family.glmm)
{
	if(is.character(family.glmm))
		family.glmm<-get(family.glmm,mode="function",envir=parent.frame())
	if(is.function(family.glmm))
		family.glmm<-family.glmm()
	#right way to check the class if there might be more than one class (maybe bc of hierarchies)
	if(!inherits(family.glmm,"glmm.family")) 
		stop(" 'family' not recognized") 
	
	family.glmm	
}

poisson.glmm <-function()
{
	family.glmm<- "poisson.glmm"
	cum <- function(eta) exp(eta)
	cp <- function(eta) exp(eta)
	cpp<-function(eta) exp(eta)
	checkData<-function(x) {
		bads<-sum(x<0)
		if(bads>0) stop("data must be nonnegative integers.")
		return(NULL)	
		}
	out<-list(family.glmm=family.glmm, cum=cum,cp=cp,cpp=cpp, checkData=checkData)
	class(out)<-"glmm.family"
	return(out)
}


bernoulli.glmm <-function()
{
	family.glmm<- "bernoulli.glmm"
	cum <- function(eta){ 
		if(eta>0) eta+log1p(exp(-eta))
		else log1p(exp(eta))}
	cum<-Vectorize(cum)
	cp <- function(eta) {1/(1+exp(-eta))}
	cpp<-function(eta) {(1/(1+exp(-eta)))*(1/(1+exp(eta)))}
	checkData<-function(x) {
		bads<-sum(x!=1&x!=0)
		if(bads>0) stop("response must be 0 and 1.")
		return(NULL)	
		}
	out<-list(family.glmm=family.glmm, cum=cum, cp=cp, cpp=cpp, checkData=checkData)
	class(out)<-"glmm.family"
	return(out)
}
