#include "myheader.h"
/*problem with gamma and squaretop causing issue with numsum: overflow/underflow*/

/*
 gamma: scalar
 thing: scalar
 squaretop: vector length m
 numsum: vector of length (nbeta + T)^2 = npar^2
 y: vector of length n
 Umat: myq by m matrix.  one COLUMN of Umat used at a time. Umat must be the R mat transposed
 myq: total number of random effects (q). scalar. 120 for salamanders.
 m: Monte Carlo sample size. scalar.
scalar.
 x: n by nbeta matrix
 n: sample size. scalar. 
 nbeta: scalar equal to the dim of beta (number of fixed effects)
 beta: vector length nbeta
 z: n by myq matrix
 Dinvfornu: square, symmetric matrix. myq x myq. D is variance matrix.
 logdetDinvfornu: scalar equal to the log determinant for D's inverse
 family_glmm:  int scalar
 Dstarinv: square, symmetric matrix. myq x myq. Dstar is variance matrix from PQL.
 logdetDstarinv: scalar equal to the log determinant for Dstar's inverse
 ustar: vector containing PQL predictions for random effects. length = myq.
 Sigmuhinv: variance matrix for normal used in importance sampling dist. myq x myq.
 logdetSigmuhinv: scalar. 
 pee: vector of length nps. vector of sampling proportions (usually 1/3, 1/3, 1/3)
 nps: scalar
 T: scalar equal to the number of variance components
 nrandom: vector of length T. 
 meow: vector of ints. length T+1. tells us which random effects go with each variance.
 nu: vector of length T. equals the variance components.
 zeta: int equal to the df for the t distrib. (p 13 of design doc)
 tconst: scalar. (p13 of design doc, equation 59)
 ntrials: is a vec of ints with length equal to n
 lfuval: scalar
 lfyuval: scalar
 wts: weights for the observations (to calculate weighted likelihood).  vector length n.

 
*/
void mcsec(double *gamma, double *thing, double *squaretop, double *numsum, double *y, double *Umat, int *myq, int *m, double *x, int *n, int *nbeta, double *beta, double *z, double *Dinvfornu, double *logdetDinvfornu, int *family_glmm, double *Dstarinv, double *logdetDstarinv, double *ustar, double *Sigmuhinv, double *logdetSigmuhinv, double *pee, int *nps, int *T, int *nrandom, int *meow, double *nu, int *zeta, double *tconst,  int *ntrials, double *lfuval, double *lfyuval, double *wts)
{
	double *Uk = Calloc(*myq, double);
	int Uindex = 0;

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
	/*double lfuval = 1.1;
	double lfyuval = 1.1;*/
	double lfutwid = 1.1;
	double *b = Calloc(*m,double);
	double a = 0.0;

	for(int k = 0; k < *m; k++){
		/*start by getting Uk  */
		for(int i = 0; i < *myq; i++){
			Uk[i] = Umat[Uindex];
			Uindex++;
		}

		/* calculate eta for this value of Uk
		first calculate ZUk*/
		matvecmult(z, Uk, n, myq, zu);

		/* then add xbeta+zu to get current value of eta */
		addvec(xbeta, zu, n, eta);

		/*log f_theta(u_k) goes into lfuval*/
		distRandGenC(Dinvfornu, logdetDinvfornu, myq, Uk, qzeros, lfuval); /*&*/

		/* log f_theta(y|u_k) value goes into lfyuval */
		elval(y, n, nbeta, eta, family_glmm, ntrials, wts, lfyuval); /*&*/

		/* value of log f~_theta(u_k) 
		first calculate value of log f~_theta(u_k) for 3 distribs used 
		and find largest value as you go */
		/* distRandGenC(Dinvfornu,logdetDinvfornu, myq, Uk, qzeros, double1);
		first piece is dist for t(0,Dstar)*/
/*		tempmax=*lfuval;*/
/*		lfutwidpieces[0]=*lfuval; */
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

		b[k] = *lfuval+ *lfyuval-lfutwid;

		if(k==0){a = b[k];}
		if(b[k]>a){a = b[k];}

		squaretop[k] = exp(2* *lfuval + 2* *lfyuval - 2*lfutwid); /*also exp(2*b[k])?*/
        /*tryfix[k] = b[k];*/
	}

	Free(lfutwidpieces);

	/* Calculate unnormalized log weights  */
	thing[0] = 0;
	for(int i = 0; i<*m; i++){
		thing[0]+= exp(b[i]-a);
	}
	Free(b);

	/* Calculate gamma */
	gamma[0] = exp(a) * *thing / *m;
    /*gamma[0] = exp(a - log(*m) +log(*thing));*/


/* now going to do second loop, for numerator */
	int npar = *nbeta + *T;
	double *lfugradient = Calloc(*T, double);
	double *lfuhess = Calloc((*T)*(*T), double);
	double *lfyugradient = Calloc(*nbeta, double);
	double *lfyuhess = Calloc((*nbeta)*(*nbeta), double);
	double *Gpiece = Calloc(npar, double);
	double *GGT = Calloc(npar*npar, double);

	Uindex = 0;	
	int Gindex = 0;
	int intone = 1;

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

		/*now I have lfugradient and lfyugradient*/

		/*create Gpiece to hold nabla log f_theta (uk,y)*/
		Gindex=0;
		for(int i = 0; i<*nbeta; i++){
			Gpiece[Gindex] = lfyugradient[i];
			Gindex++;
		}
		for(int i = 0; i<*T; i++){
			Gpiece[Gindex] = lfugradient[i];
			Gindex++;
		}
		/*now Gpiece is ready*/

		/*calculate GGT Gpiece t(Gpiece)*/
		matmatmult(Gpiece, Gpiece, &npar, &intone, &npar, GGT);
		
		/*add on GGT * weight to numsum. this corresponds to numpieces from R*/
		for(int i = 0; i < npar*npar; i++){
			numsum[i]+= GGT[i] * squaretop[k];
		}

	} /* ends  k loop */

	Free(GGT);
	Free(Gpiece);
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


