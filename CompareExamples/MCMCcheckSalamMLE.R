library(bernor,lib.loc="bernor.Rcheck")
 library(mcmc)

 data(salam)

 ##### global variables #####

 attach(salam) # contains x, y, z, i

 nparm <- ncol(x) + length(unique(i))

 ##### seeds #####

 set.seed(42)

 ##### functions #####
#calculates value of sum CLL
 h <- function(b) {
     b <- matrix(b, ncol = ncol(y))
     result <- 0
     for (j in 1:ncol(y))
         result <- result + bernor(y[ , j], mu, b[ , j], sigma, x, z, i)$value
     return(result)
 }

#calculates gradient of sum CLL
 g <- function(b) {
     b <- matrix(b, ncol = ncol(y))
     result <- rep(0, nparm)
     for (j in 1:ncol(y))
         result <- result + bernor(y[ , j], mu, b[ , j], sigma, x, z, i,
             deriv = 1)$gradient
     return(result)
 }

##### try 1 #####

 # from Booth and Hobert (1999)
#mu<- RR, RW, WR, WW
#sigma<- sqrt(c(nuF, nuM))
mu<-c(1,.3,-2,1)
 sigma <- sqrt(c(1.5, 1.3))

 b <- rep(0, length(i) * ncol(y))
 out <- metrop(h, b, nbatch = 100, blen = 1000, outfun = g)
 print(out$time)
 print(out$accept)

#try 2
 out <- metrop(out,scale=.1)
 print(out$time)
 print(out$accept)

 print(apply(out$batch, 2, mean))
 print(apply(out$batch, 2, sd) / sqrt(out$nbatch))

#looks good enough so collect enough u's (b's here)

out<-metrop(out,blen=10^3)
print(out$accept)
 print(apply(out$batch, 2, mean))
 print(apply(out$batch, 2, sd) / sqrt(out$nbatch))

themean<- print(apply(out$batch, 2, mean))
ME<-2* print(apply(out$batch, 2, sd) / sqrt(out$nbatch))
cbind(themean-ME,themean+ME)



