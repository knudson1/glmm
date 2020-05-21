#include "myheader.h"
/*
 y: vector of length n
 Umat: myq by m matrix, one COLUMN of Umat used at a time. Umat must be the R mat transposed
 myq: total number of random effects (q). scalar. 120 for salamanders.
 m: Monte Carlo sample size. scalar.		 
 x: n by nbeta matrix
 n: sample size. scalar. 
 nbeta: scalar equal to the dim of beta (number of fixed effects)
 beta: vector of length nbeta
 z:  n by myq matrix
 Dinvfornu: square, symmetric matrix. myq x myq. D is variance matrix.
 logdetDinvfornu: scalar equal to the log determinant for D's inverse
 family_glmm:  int scalar
 Dstarinv: square, symmetric matrix. myq x myq. Dstar is variance matrix from PQL.
 logdetDstarinv: scalar equal to the log determinant for Dstar's inverse
 ustar: vector containing PQL predictions for random effects. length = myq.
 Sigmuhinv: variance matrix for normal used in importance sampling dist. myq x myq.
 logdetSigmuhinv: scalar. 
 pee: vector of sampling proportions (usually 1/3, 1/3, 1/3)
 nps: scalar equal to the length of pee (3 for now, maybe more if imp sampling distrib changes)
 T: scalar equal to the number of variance components
 nrandom: vector of length T. each entry equals number of random effects associated with that variance component. For salamanders, 60 F and 60 M so nrandom is vec of length 2 (60, 60)
 meow: vector of ints. length T+1. tells us which random effects go with each variance.
 nu: vector of length T. equals the variance components.
 zeta: int equal to the df for the t distrib. (p 13 of design doc)
 tconst: scalar. (p13 of design doc, equation 59)
 v: vector of length n. normalized importance sampling weights
 ntrials:  a vec of ints with length n
 value:  scalar equal to the value of the MC log likelihood
 gradient: vector equal to the first deriv of the MC log likelihood (length nbeta + T)
 b: vector of length m. unnormalized importance sampling weights
 wts: weights for the observations (to calculate weighted likelihood)

 */

void valgrad(double *y, double *Umat, int *myq, int *m, double *x, int *n, int *nbeta, double *beta, double *z, double *Dinvfornu, double *logdetDinvfornu, int *family_glmm, double *Dstarinv, double *logdetDstarinv, double *ustar, double *Sigmuhinv, double *logdetSigmuhinv, double *pee, int *nps, int *T, int *nrandom, int *meow, double *nu, int *zeta, double *tconst, double *v, int *ntrials, double *value, double *gradient, double *b, double *wts)
{
    double *Uk = Calloc(*myq, double);
    
    double db1 = 0.0;
    double *double1 = &db1; /*temp to hold info that will be put into lfuval*/
    
    /* Calculate xbeta, needed to calculate eta for each Uk=U[k,] in R notation */
    double *xbeta = Calloc(*n,double);
    matvecmult(x,beta,n,nbeta,xbeta);
    double *zu = Calloc(*n,double);
    double *eta = Calloc(*n,double);
    double *qzeros = Calloc(*myq,double);
    double tempmax = 1.0;
    double *lfutwidpieces = Calloc(*nps,double);
    double diffs = 0.0;
    double lfuval = 1.1;
    double lfyuval = 1.1;
    double lfutwid = 1.1;
    double a = 0.0;
    
    for(int k = 0, Uindex = 0; k < *m; k++){
        /*start by getting Uk */
        for(int i = 0; i < *myq; i++, Uindex++){
            Uk[i] = Umat[Uindex];
        }
        
        /* calculate eta for this value of Uk
         first calculate ZUk*/
        matvecmult(z, Uk, n, myq, zu);
        
        /* then add xbeta+zu to get current value of eta */
        addvec(xbeta, zu, n, eta);
        
        /*log f_theta(u_k) goes into lfuval*/
        distRandGenC(Dinvfornu, logdetDinvfornu, myq, Uk, qzeros, &lfuval);
        
        /* log f_theta(y|u_k) value goes into lfyuval */
        elval(y, n, nbeta, eta, family_glmm, ntrials, wts, &lfyuval);
        
        /* value of log f~_theta(u_k)
         first calculate value of log f~_theta(u_k) for 3 distribs used
         and find largest value as you go */
        /* distRandGenC(Dinvfornu,logdetDinvfornu, myq, Uk, qzeros, double1);
         first piece is dist for t(0,Dstar)*/
        /*        tempmax=*lfuval;*/
        /*        lfutwidpieces[0]=*lfuval; */
        tdist(Dstarinv, myq, Uk, zeta, tconst, double1);
        lfutwidpieces[0] = *double1;
        tempmax = *double1;
        
        /*second piece is for N(ustar,Dstar)  */
        distRandGenC(Dstarinv,logdetDstarinv, myq, Uk, ustar, double1);
        lfutwidpieces[1] = *double1;
        if(*double1 > tempmax){tempmax = *double1;}
        
        /*third piece is for N(ustar,complicated var)  */
        distRandGenC(Sigmuhinv,logdetSigmuhinv, myq, Uk, ustar, double1);
        lfutwidpieces[2] = *double1;
        if(*double1>tempmax){tempmax = *double1;}
        
        lfutwid=0.0;
        
        /* calculate diffs= lfutwidpieces-tempmax */
        /* and lfutwid = tempmax + log(sum(pee[i]*exp(diffs))) */
        for(int i = 0; i<*nps; i++){
            diffs = lfutwidpieces[i]- tempmax;
            lfutwid+= pee[i]*exp(diffs); /*sum(pee[i]*exp(diffs))*/
        }
        lfutwid = log(lfutwid)+tempmax; /* finishes lfutwid calc */
        
        b[k] = lfuval+lfyuval-lfutwid;
        
        if(k==0){a = b[k];}
        if(b[k]>a){a = b[k];}
    }
    
    Free(lfutwidpieces);
    
    /* Calculate weights v[k] */
    double *tops = Calloc(*m, double);
    for(int i = 0; i<*m; i++){
        tops[i] = exp(b[i]-a);
    }
    
    double bottom = 0.0;
    for(int i = 0; i<*m; i++){
        bottom+= tops[i];
    }
    /* Calculate tops/bottom */
    for(int i = 0; i<*m; i++){
        v[i] = tops[i]/(bottom);
    }
    Free(tops);
    
    /* Calculate value */
    *value = a-log(*m)+log(bottom);
    
    
    /* done with value! */
    /* now going to do second loop, which calculates grad and hess */

    double *lfugradient = Calloc(*T, double);
    double *lfuhess = Calloc((*T)*(*T), double);
    double *lfyugradient = Calloc(*nbeta, double);
    double *lfyuhess = Calloc((*nbeta)*(*nbeta), double);
    
    
    int Gindex = 0;

    
    for(int k = 0, Uindex = 0; k < *m; k++){
        /*start by getting Uk  */
        for(int i = 0; i<*myq; i++, Uindex++){
            Uk[i] = Umat[Uindex];
        }
        
        /* calculate eta for this value of Uk
         first calculate ZUk*/
        matvecmult(z, Uk, n, myq, zu);
        
        /* then add xbeta+zu to get current value of eta */
        addvec(xbeta, zu, n, eta);
        
        /* calculate lfu gradient and hessian */
        distRand3C(nu, qzeros, T, nrandom, meow, Uk, lfugradient, lfuhess);
        
        /* calculate gradient and hessian log f_theta(y|u_k) */
        elGH(y, x, n, nbeta, eta, family_glmm, ntrials, wts, lfyugradient, lfyuhess);
        
        /* calculate gradient */
        Gindex=0;
        
        for(int i = 0; i<*nbeta; i++){
            gradient[Gindex]+= lfyugradient[i]*v[k];
            Gindex++;
        }
        for(int i = 0; i<*T; i++){
            gradient[Gindex]+= lfugradient[i]*v[k];
            Gindex++;
        }
        
    } /* ends FIRST k loop */
    
    Free(Uk);
    Free(xbeta);
    Free(zu);
    Free(eta);
    Free(lfugradient);
    Free(lfuhess);
    Free(lfyugradient);
    Free(lfyuhess);
    Free(qzeros);
    
}


