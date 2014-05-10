getFamily <-
function(family.glmm)
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

poisson.glmm <-
function()
{
	family.glmm<- "poisson.glmm"
	cum <- function(eta) exp(eta)
	cp <- function(eta) exp(eta)
	cpp<-function(eta) exp(eta)
	checkData<-function(x) {
		bads<-sum(x<0)
		if(bads>0) stop("data must be nonnegative.")
		return(NULL)	
		}
	out<-list(family.glmm=family.glmm, cum=cum,cp=cp,cpp=cpp, checkData=checkData)
	class(out)<-"glmm.family"
	return(out)
}

bernoulli.glmm <-
function()
{
	family.glmm<- "bernoulli.glmm"
	cum <- function(eta) log(1+exp(eta))
	cp <- function(eta) exp(eta)/(1+exp(eta))
	cpp<-function(eta) exp(eta)/(1+exp(eta))-exp(2*eta)*(1+exp(eta))^(-2)
	
	checkData<-function(x) {
		bads<-sum(x>1) +sum(x<0)
		if(bads>0) stop("data must be between 0 and 1.")
		return(NULL)	
		}
	out<-list(family.glmm=family.glmm, cum=cum, cp=cp, cpp=cpp, checkData=checkData)
	class(out)<-"glmm.family"
	return(out)
}
