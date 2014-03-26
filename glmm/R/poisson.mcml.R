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
