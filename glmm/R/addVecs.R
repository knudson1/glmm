addVecs <-
function(vecs){
	 if (! is.list(vecs))
        vecs <- list(vecs)

	dvecs <-length(vecs[[1]])
	out<-rep(0,dvecs)
	
	
	for(d in 1:dvecs){
		thing<-lapply(vecs,"[[",d)
		thing<-unlist(thing)
		out[d]<-sum(thing)
	}
	out
}
