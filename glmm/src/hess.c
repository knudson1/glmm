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
 gradient: vector equal to the first deriv of the MC log likelihood (length nbeta + T)
 hessian: square matrix with nrows = (nbeta+T) 
 b: vector of length m. unnormalized importance sampling weights
 length: 
 q: 
 wts: weights for the observations (to calculate weighted likelihood)

 */
// value can go?

// Charlie: shut warnings about unused parameters only way CRAN allows
#if defined(__GNUC__) || defined(__clang__)
void hess(double *y, 
double *Umat, 
int *myq, 
int *m, 
double *x, 
int *n, 
int *nbeta, 
double *beta, 
double *z, 
double *Dinvfornu __attribute__ ((unused)), 
double *logdetDinvfornu __attribute__ ((unused)), 
int *family_glmm, 
double *Dstarinv __attribute__ ((unused)), 
double *logdetDstarinv __attribute__ ((unused)),
 double *ustar __attribute__ ((unused)), 
double *Sigmuhinv __attribute__ ((unused)), 
double *logdetSigmuhinv __attribute__ ((unused)), 
double *pee __attribute__ ((unused)), 
int *nps __attribute__ ((unused)), 
int *T, 
int *nrandom, 
int *meow, 
double *nu, 
int *zeta __attribute__ ((unused)), 
double *tconst __attribute__ ((unused)), 
double *v, 
int *ntrials, 
double *gradient, 
double *hessian, 
double *b, 
int *length, 
double *q, 
double *wts)
#else
void hess(double *y, double *Umat, int *myq, int *m, double *x, int *n, int *nbeta, double *beta, double *z, double *Dinvfornu, double *logdetDinvfornu, int *family_glmm, double *Dstarinv, double *logdetDstarinv, double *ustar, double *Sigmuhinv, double *logdetSigmuhinv, double *pee, int *nps, int *T, int *nrandom, int *meow, double *nu, int *zeta, double *tconst, double *v, int *ntrials, double *gradient, double *hessian, double *b, int *length, double *q, double *wts)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    double *Uk = R_Calloc(*myq, double);
    // clang static analyzer complains about the following
    // put the int on the assignment below that ignores this one
    // int Uindex = 0;
    
    /* Calculate xbeta, needed to calculate eta for each Uk=U[k,] in R notation */
    double *xbeta = R_Calloc(*n,double);
    matvecmult(x,beta,n,nbeta,xbeta);
    double *zu = R_Calloc(*n,double);
    double *eta = R_Calloc(*n,double);
    double *qzeros = R_Calloc(*myq,double);
    double a = 0.0;
    
    for(int k = 0; k < *m; k++){
        
        if(k==0){a = b[k];}
        if(b[k]>a){a = b[k];}
    }
    
    
    /* Calculate weights v[k] */
    double *tops = R_Calloc(*m, double);
    for(int i = 0; i<*m; i++){
        tops[i] = exp(q[i]-a);
    }
    
    double *bs = R_Calloc(*length, double);
    for(int i = 0; i<*length; i++){
        bs[i] = exp(b[i]-a);
    }
    
     double bottom = 0.0;
     for(int i = 0; i<*length; i++){
        bottom+= bs[i];
    }
    /* Calculate tops/bottom */
    for(int i = 0; i<*m; i++){
        v[i] = tops[i]/(bottom);
    }
    R_Free(tops);
    
    /* done with value! */
    /* now going to do second loop, which calculates grad and hess */
    int npar = *nbeta + *T;
    double *lfugradient = R_Calloc(*T, double);
    double *lfuhess = R_Calloc((*T)*(*T), double);
    double *lfyugradient = R_Calloc(*nbeta, double);
    double *lfyuhess = R_Calloc((*nbeta)*(*nbeta), double);
    
    double *pandabit1 = R_Calloc(npar, double);
    double *pandabit2 = R_Calloc(npar, double);
    double *panda = R_Calloc(npar*npar, double);
    double *pandatemp = R_Calloc(npar*npar, double);
    
    int Uindex = 0;
    int Gindex = 0;
    int intone = 1;
    
    double *lobster = R_Calloc(npar*npar, double);
    int lfyuindex = 0, lfuindex = 0, matindex = 0;
    
    
    
    /*    reset counters to 0*/
    Uindex = 0;
    // Gindex = 0;     // not needed; reset in loop
    // lfyuindex = 0;  // not needed; reset in loop
    // lfuindex = 0;   // not needed; reset in loop
    // matindex = 0;   // not needed; reset in loop
    
    /* begins SECOND k loop */
    for(int k = 0; k<*m; k++){
        /*start by getting Uk  */
        for(int i = 0; i<*myq; i++){
            Uk[i] = Umat[Uindex];
            Uindex++;
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
        
        /* Calculate hessian. 2 parts: panda, lobster.
         Panda is \sum_{k=1}^m pandabit1 times pandabit2 where
         pandabit1 is (\nabla log f_\theta (uk,y)-\nabla objfun)
         and pandabit2 is (\nabla log f_\theta (uk,y)-\nabla objfun)' v(uk,y)
         Lobster is \sum_{k=1}^m \nabla^2 log f_\theta(uk,y) v(uk,y)
         */
        Gindex = 0;
        for(int i = 0; i<*nbeta; i++){
            pandabit1[Gindex] = lfyugradient[i] - gradient[Gindex];
            pandabit2[Gindex] = (lfyugradient[i] - gradient[Gindex])*v[k];
            Gindex++;
        }
        for(int i = 0; i<*T; i++){
            pandabit1[Gindex] = lfugradient[i] - gradient[Gindex];
            pandabit2[Gindex] = (lfugradient[i] - gradient[Gindex])*v[k];
            Gindex++;
        }
        /* now have pandabits so can calc what we'll add on this iteration*/
        matmatmult(pandabit1, pandabit2, &npar, &intone, &npar, pandatemp);
        
        /* now do panda+=pandatemp */
        for(int i = 0; i<(npar*npar); i++){
            panda[i]+= pandatemp[i];
        }
        
        /* Calculate lobster */
        lfyuindex = 0; /* every k-iteration want to start at */
        lfuindex = 0; /* beginning of array for lfyuhess, lfuhess, lobster */
        matindex = 0;
        for(int i = 0; i<*nbeta; i++){ /* for first nbeta columns of lobster */
            for(int j = 0; j<*nbeta; j++){
                lobster[matindex]+= lfyuhess[lfyuindex]*v[k];
                matindex++;
                lfyuindex++;
            }
            for(int j = 0; j<*T; j++){
                matindex++; /* add 0 to lobster */
            }
        }
        for(int i = 0; i<*T; i++){ /* for last T cols of lobster */
            for(int j = 0; j<*nbeta; j++){
                matindex++; /* add 0 to lobster */
            }
            for(int j = 0; j<*T; j++){
                lobster[matindex]+= lfuhess[lfuindex]*v[k];
                matindex++;
                lfuindex++;
            } /* ends j loop */
        }/* ends i loop */
    } /* ends SECOND k loop */
    
    R_Free(pandatemp);
    R_Free(pandabit1);
    R_Free(pandabit2);
    R_Free(Uk);
    R_Free(xbeta);
    R_Free(zu);
    R_Free(eta);
    R_Free(lfugradient);
    R_Free(lfuhess);
    R_Free(lfyugradient);
    R_Free(lfyuhess);
    R_Free(qzeros);
    
    /* finally, hessian = lobster + panda  */
    for(int i = 0; i<(npar*npar); i++){
        hessian[i] = lobster[i] + panda[i] ;
    }
    
    R_Free(panda);
    R_Free(lobster);

    // caught by ../../devel/calloc-match.R
    R_Free(bs);

}


