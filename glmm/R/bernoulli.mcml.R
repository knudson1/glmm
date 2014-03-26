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
