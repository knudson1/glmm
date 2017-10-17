library(glmm,lib.loc="glmm.Rcheck")
set.seed(1234)
data(salamander)
m<-10^4
sal<-glmm(Mate~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"), data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE)
summary(sal)

