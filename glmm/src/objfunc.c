#include "myheader.h"
/*x is n by nbeta matrix
beta has length nbeta
Umat is myq by m matrix, one COLUMN of Umat used at a time. Umat must be the R mat transposed 
z is n by myq matrix
pee is the vector of sampling proportions (usually 1/3, 1/3, 1/3)
nps is the length of pee (3 for now, maybe more if imp sampling distrib changes)
*/
void objfunc(double *y, double *Umat, int *myq, int *m, double *x, int *n, int *nbeta, double *beta, double *z, double *Dinvfornu, double *logdetDinvfornu, int *family_glmm, double *Dstarinv, double *logdetDstarinv, double *ustar, double *Sigmuhinv, double *logdetSigmuhinv, double *pee, int *nps, int *T, int *nrandom, int *meow, double *nu, int *zeta, double *tconst, double *v, double *value, double *gradient, double *hessian)
{
	double *Uk=Calloc(myq,double);
	int Uindex=0;

	double db1=0.0;
	double *double1=&db1; /*temp to hold info that will be put into lfuval*/

	/* Calculate xbeta, needed to calculate eta for each Uk=U[k,] in R notation */
	double *xbeta=Calloc(*n,double);
	matvecmult(x,beta,n,nbeta,xbeta);
	double *zu=Calloc(*n,double);
	double *eta=Calloc(*n,double);
	double *qzeros=Calloc(*myq,double);
	double tempmax=1.0;
	double *lfutwidpieces=Calloc(*nps,double);
	double diffs=0.0;
	double *lfuval=Calloc(1,double);
	double *lfyuval=Calloc(1,double);
	double *lfutwid=Calloc(1,double);
	double *b=Calloc(*m,double);
	double *a=Calloc(1,double);

	for(int k=0;k<*m;k++){
		/*start by getting Uk  */
		for(int i=0;i<*myq;i++){
			Uk[i]=Umat[Uindex];
			Uindex++;
		}

		/* calculate eta for this value of Uk
		first calculate ZUk*/
		matvecmult(z,Uk,n,myq,zu);

		/* then add xbeta+zu to get current value of eta */
		addvec(xbeta,zu,n,eta);

		/*log f_theta(u_k) goes into lfuval*/
		distRandGenC(Dinvfornu,logdetDinvfornu,myq,Uk,qzeros,lfuval); 

		/* log f_theta(y|u_k) value goes into lfyuval */
		elval(y,n,nbeta,eta,family_glmm,lfyuval);

		/* value of log f~_theta(u_k) 
		first calculate value of log f~_theta(u_k) for 3 distribs used 
		and find largest value as you go */
		/* distRandGenC(Dinvfornu,logdetDinvfornu, myq, Uk, qzeros, double1);
		first piece is dist for t(0,Dstar)*/
/*		tempmax=*lfuval;*/
/*		lfutwidpieces[0]=*lfuval; */
		tdist(Dstarinv, myq, Uk, zeta, tconst, double1);
		lfutwidpieces[0]=*double1;
		tempmax=*double1;
		
		/*second piece is for N(ustar,Dstar)  */
		distRandGenC(Dstarinv,logdetDstarinv, myq, Uk, ustar, double1);
		lfutwidpieces[1]=*double1;
		if(*double1>tempmax){tempmax=*double1;}

		/*third piece is for N(ustar,complicated var)  */
		distRandGenC(Sigmuhinv,logdetSigmuhinv, myq, Uk, ustar, double1);
		lfutwidpieces[2]=*double1;
		if(*double1>tempmax){tempmax=*double1;}

		*lfutwid=0.0;

		/* calculate diffs= lfutwidpieces-tempmax */
		/* and lfutwid = tempmax + log(sum(pee[i]*exp(diffs))) */
		for(int i=0;i<*nps;i++){
			diffs=lfutwidpieces[i]- tempmax;
			*lfutwid+=pee[i]*exp(diffs); /*sum(pee[i]*exp(diffs))*/
		}
	*lfutwid=log(*lfutwid)+tempmax; /* finishes lfutwid calc */

	b[k]=*lfuval+*lfyuval-*lfutwid;

	if(k==0){*a=b[k];}
	if(b[k]>*a){*a=b[k];}
	}

	Free(lfutwidpieces);
	Free(lfuval);
	Free(lfyuval);
	Free(lfutwid);

	/* Calculate weights v[k] */
	double *tops=Calloc(*m,double);
	for(int i=0;i<*m;i++){
		tops[i]=exp(b[i]-*a);
	}
	Free(b);

	double bottom=0.0;
	for(int i=0;i<*m;i++){
		bottom+=tops[i];
	}
	/* Calculate tops/bottom */
	for(int i=0;i<*m;i++){
		v[i]=tops[i]/bottom;
	}
	Free(tops);

	/* Calculate value */
	*value=*a-log(*m)+log(bottom);
	Free(a);

/* done with value! */
/* now going to do second loop, which calculates grad and hess */
	int npar=*nbeta+*T;
	double *lfugradient=Calloc(*T,double);
	double *lfuhess=Calloc((*T)*(*T),double);
	double *lfyugradient=Calloc(*nbeta,double);
	double *lfyuhess=Calloc((*nbeta)*(*nbeta),double);

	double *pandabit1=Calloc(npar,double);
	double *pandabit2=Calloc(npar,double);
	double *panda=Calloc(npar*npar,double);
	double *pandatemp=Calloc(npar*npar,double);

	Uindex=0;	
	int Gindex=0;
	int intone=1;

	double *lobster=Calloc(npar*npar,double);
	int lfyuindex=0,lfuindex=0,matindex=0;

	for(int k=0;k<*m;k++){
		/*start by getting Uk  */
		for(int i=0;i<*myq;i++){
			Uk[i]=Umat[Uindex];
			Uindex++;
		}

		/* calculate eta for this value of Uk
		first calculate ZUk*/
		matvecmult(z,Uk,n,myq,zu);

		/* then add xbeta+zu to get current value of eta */
		addvec(xbeta,zu,n,eta);

		/* calculate lfu gradient and hessian */
		distRand3C(nu, qzeros, T, nrandom, meow, Uk, lfugradient, lfuhess);

		/* calculate gradient and hessian log f_theta(y|u_k) */
		elGH(y,x,n,nbeta,eta,family_glmm,lfyugradient,lfyuhess);

		/* calculate gradient */
		Gindex=0;

		for(int i=0;i<*nbeta;i++){
			gradient[Gindex]+=lfyugradient[i]*v[k];
			Gindex++;
		}
		for(int i=0;i<*T;i++){
			gradient[Gindex]+=lfugradient[i]*v[k];
			Gindex++;
		}

	/* Calculate hessian. 3 parts: panda, lobster, and matrix G t(G). 
	Panda is \sum_{k=1}^m (\nabla log f_\theta (uk,y))(\nabla log f_\theta (uk,y))' v(uk,y)
	Lobster is \sum_{k=1}^m \nabla^2 log f_\theta(uk,y) v(uk,y)  
	panda and lobster are calculated in a loop, but G t(G) will be calc outside of k-loop.
		First calculate panda. For each k, add on pandabit1 %*% t(pandabit2) 
		where pandabit2=pandabit1* v[k]. */
		Gindex=0;
		for(int i=0;i<*nbeta;i++){
			pandabit1[Gindex]=lfyugradient[i];
			pandabit2[Gindex]=lfyugradient[i]*v[k];
			Gindex++;
		}
		for(int i=0;i<*T;i++){
			pandabit1[Gindex]=lfugradient[i];
			pandabit2[Gindex]=lfugradient[i]*v[k];
			Gindex++;
		}
		/* now have pandabits so can calc what we'll add on this iteration*/
		matmatmult(pandabit1,pandabit2,&npar,&intone,&npar,pandatemp);
		
		/* now do panda+=pandatemp */
		for(int i=0;i<(npar*npar);i++){
			panda[i]+=pandatemp[i];
		}

		/* Calculate lobster */
		lfyuindex=0; /* every k-iteration want to start at */
		lfuindex=0; /* beginning of array for lfyuhess, lfuhess, lobster */
		matindex=0;
		for(int i=0;i<*nbeta;i++){ /* for first nbeta columns of lobster */
			for(int j=0;j<*nbeta;j++){
				lobster[matindex]+=lfyuhess[lfyuindex]*v[k];
				matindex++;
				lfyuindex++;
			}
			for(int j=0;j<*T;j++){
				matindex++; /* add 0 to lobster */
			}
		}
		for(int i=0;i<*T;i++){ /* for last T cols of lobster */
			for(int j=0;j<*nbeta;j++){
				matindex++; /* add 0 to lobster */
			}
			for(int j=0;j<*T;j++){
				lobster[matindex]+=lfuhess[lfuindex]*v[k];
				matindex++;
				lfuindex++;
			}
		}
	}

	Free(pandatemp);
	Free(pandabit1);
	Free(pandabit2);
	Free(Uk);
	Free(xbeta);
	Free(zu);
	Free(eta);
	Free(lfugradient);
	Free(lfuhess);
	Free(lfyugradient);
	Free(lfyuhess);
	Free(qzeros);

	double *GtG=Calloc((npar*npar),double);
	matmatmult(gradient,gradient,&npar,&intone,&npar,GtG);

	/* finally, hessian = lobster+ panda - gradient t(gradient) */
	for(int i=0;i<(npar*npar);i++){
		hessian[i]=lobster[i]+panda[i]-GtG[i];
	}

	Free(panda);
	Free(GtG);
	Free(lobster);

}


