addMats <-
function(matList){
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
