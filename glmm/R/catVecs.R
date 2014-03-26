catVecs <-
function(vecs){
	 if (! is.list(vecs))
        vecs <- list(vecs)

	lengths <-lapply(vecs,length)
	lengths<-unlist(lengths)
	out<-rep(0,sum(lengths))
	starthere<-1
	
	for(d in 1:length(vecs)){
		endhere<-starthere+lengths[d]-1
		out[c(starthere:endhere)]<-unlist(vecs[[d]])
		starthere<-starthere+lengths[d]
	}
	out
}
