print.summary.mcla <-
function(x,...){
    summ<-x	
    stopifnot(inherits(summ, "summary.mcla"))

    cat("\nCall:\n", paste(deparse(summ$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
   
cat("Fixed Effects:")
   cat("\n")

	printCoefmat(summ$coefmat,digits=3)
   cat("\n")

   cat("\n")
cat("Variance Components for Random Effects (P-values are one-tailed):")
   cat("\n")

	printCoefmat(summ$nucoefmat,digits=3)
   cat("\n")

}
