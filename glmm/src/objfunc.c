#include "myheader.h"
/*x is n by nbeta matrix
beta has length nbeta
Umat is myq by m matrix, one COLUMN of Umat used at a time. Umat must be the R mat transposed 
z is n by myq matrix
pee is the vector of sampling proportions (usually 1/3, 1/3, 1/3)
nps is the length of pee (3 for now, maybe more if imp sampling distrib changes)
*/
void objfunc(double *y, double *Umat, int *myq, int *m, double *x, int *n, int *nbeta, double *beta, double *z, double *Dinvfornu, double *logdetDinvfornu, int *family_glmm, double *Dstarinv, double *logdetDstarinv, double *ustar, double *Sigmuhinv, double *logdetSigmuhinv, double *pee, int *nps, int *T, int *nrandom, int *meow, double *nu, double *v, double *value, double *G, double *hessian)
{
	double *Uk=Calloc(myq,double);
	int k=0,i=0,Uindex=0;

	double db1=0.0,*double1=&db1; /*temp to hold info that will be put into lfuval*/

	/* calculate xbeta, used to calculate eta for each Uk=U[k,] in R notation */
	double *xbeta=Calloc(*n,double);
	matvecmult(x,beta,n,nbeta,xbeta);
	double *zu=Calloc(*n,double);
	double *eta=Calloc(*n,double);
	double *mzeros=Calloc(*m,double);
	memset(mzeros,0,*m);
	double tempmax=1.0;
	double *lfutwidpieces=Calloc(*nps,double);
	double diffs=0.0;
	double *lfuval=Calloc(*m,double);
	double *lfyuval=Calloc(*m,double);
	double *lfutwid=Calloc(*m,double);
	double *b=Calloc(*m,double);
	double *a=Calloc(1,double);

	for(k=0;k<*m;k++){
		/*start by getting Uk  */
		for(i=0;i<*myq;i++){
			Uk[i]=Umat[Uindex];
			Uindex++;
		}

		/* calculate eta for this value of Uk
		first calculate ZUk*/
		matvecmult(z,Uk,n,myq,zu);

		/* then add xbeta+zu to get current value of eta */
		addvec(xbeta,zu,n,eta);

		/*log f_theta(u_k)*/
		distRandGenC(Dinvfornu,logdetDinvfornu,myq,Uk,mzeros,double1); 
		*(lfuval+k)=*double1; 

		/* log f_theta(y|u_k) value */
		elval(y,x,n,nbeta,eta,family_glmm,double1);
		*(lfyuval+k)=*double1;

		/* value of log f~_theta(u_k) 
		first calculate value of log f~_theta(u_k) for 3 distribs used 
		and find largest value as you go */
		distRandGenC(Dstarinv,logdetDstarinv, myq, Uk, mzeros, double1);
		tempmax=*double1;
		lfutwidpieces[0]=*double1;
		
		distRandGenC(Dstarinv,logdetDstarinv, myq, Uk, ustar, double1);
		lfutwidpieces[1]=*double1;
		if(*double1>tempmax){tempmax=*double1;}

		distRandGenC(Sigmuhinv,logdetSigmuhinv, myq, Uk, ustar, double1);
		lfutwidpieces[2]=*double1;
		if(*double1>tempmax){tempmax=*double1;}

		lfutwid[k]=0.0;

		/* calculate diffs= lfutwidpieces-tempmax */
		/* and lfutwid = tempmax + log(sum(pee[i]*exp(diffs))) */
		for(i=0;i<*nps;i++){
			diffs=lfutwidpieces[i]- tempmax;
			lfutwid[k]+=*(pee+i)*exp(diffs); /*sum(pee[i]*exp(diffs))*/
		}
	lfutwid[k]=log(lfutwid[k])+tempmax; /* finishes lfutwid calc */

	b[k]=*(lfuval+k)+*(lfyuval+k)-*(lfutwid+k);

	if(k==0){*a=b[k];}
	if(b[k]>*a){*a=b[k];}
	}


	Free(lfutwidpieces);
	Free(lfuval);
	Free(lfyuval);
	Free(lfutwid);

	/* Calculate weights v[k] */
	double *tops=Calloc(*m,double);
	for(i=0;i<*m;i++){
		tops[i]=exp(b[i]-*a);
	}
	Free(b);

	double bottom=0.0;
	for(i=0;i<*m;i++){
		bottom+=tops[i];
	}
	/* calc tops/bottom */
	for(i=0;i<*m;i++){
		v[i]=tops[i]/bottom;
	}
	Free(tops);

	/* calculate value */
	*value=*a-log(*m)+log(bottom);
	Free(a);

/* done with value! */
/* now going to do second loop, which calcs grad and hess */
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
	int thing1=1,*ione=&thing1;

	for(k=0;k<*m;k++){
		/*start by getting Uk  */
		for(i=0;i<*myq;i++){
			Uk[i]=Umat[Uindex];
			Uindex++;
		}

		/* calculate eta for this value of Uk
		first calculate ZUk*/
		matvecmult(z,Uk,n,myq,zu);

		/* then add xbeta+zu to get current value of eta */
		addvec(xbeta,zu,n,eta);

		/* calculate lfu gradient and hessian */
		distRand3C(nu, mzeros, T, nrandom, meow, Uk, lfugradient, lfuhess);

		/* calculate gradient and hessian log f_theta(y|u_k) */
		elGH(y,x,n,nbeta,eta,family_glmm,lfyugradient,lfyuhess);

		/* calculate G */
		Gindex=0;

		for(i=0;i<*nbeta;i++){
			G[Gindex]+=lfyugradient[i]*v[k];
			Gindex++;
		}
		for(i=0;i<*T;i++){
			G[Gindex]+=lfugradient[i]*v[k];
			Gindex++;
		}

	/* Now moving on to HESSIAN. 3 parts: panda, lobster, and matrix G t(G) 
	need loops for panda and lobster, but G t(G) will be calc at end.
		First do panda. for each k, add on pandabit1 %*% t(pandabit2) 
		where pandabit2=pandabit1* v[k]. */
		Gindex=0;
		for(i=0;i<*nbeta;i++){
			pandabit1[Gindex]=lfyugradient[i];
			pandabit2[Gindex]=lfyugradient[i]*v[k];
			Gindex++;
		}
		for(i=0;i<*T;i++){
			pandabit1[Gindex]=lfugradient[i];
			pandabit2[Gindex]=lfugradient[i]*v[k];
			Gindex++;
		}
		/* now have pandabits so can calc what we'll add on this iteration*/
		matmatmult(pandabit1,pandabit2,&npar,ione,&npar,pandatemp);
		
		/* now do panda+=pandatemp */
		for(i=0;i<(npar*npar);i++){
			panda[i]+=pandatemp[i];
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
	Free(mzeros);

	/* finally, hessian = lobster+ panda + G t(G) */
	for(i=0;i<(npar*npar);i++){
		hessian[i]=panda[i];
	}

	Free(panda);

}


