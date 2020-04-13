# This is the examples section from logLik.Rd

library(glmm, lib.loc = "../glmm.Rcheck")
data(BoothHobert)
set.seed(1234)

mod <- glmm(y~0+x1, list(y~0+z1), varcomps.names=c("z1"), 
    data=BoothHobert, family.glmm=bernoulli.glmm, m=100)
logLik(mod)
