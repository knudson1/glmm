library(glmm)
set.seed(1234)

 nrow <- 5
 ncol <- 3

 a <- matrix(rnorm(nrow * ncol), nrow = nrow)
 b <- rnorm(ncol)

 foo1 <- as.numeric(a %*% b)

 mout3 <- .C("matvecmult", a = as.double(a), b = as.double(b),
     nrow = as.integer(nrow), ncol = as.integer(ncol), result = double(nrow))
 identical(foo1, mout3$result)
 all.equal(foo1, mout3$result)

mout3









