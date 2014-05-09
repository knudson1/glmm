getFamily <-
function(family.mcml)
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

poisson.mcml <-
function()
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

bernoulli.mcml <-
function()
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
