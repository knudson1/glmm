#include "myheader.h"
/*x is n by nbeta matrix
beta has length nbeta
Umat is myq by m matrix, one COLUMN of Umat used at a time. Umat must be the R mat transposed 
z is n by myq matrix
pee is the vector of sampling proportions (usually 1/3, 1/3, 1/3)
nps is the length of pee (3 for now, maybe more if imp sampling distrib changes)
*/
void objfunc(double *y, double *Umat, int *myq, int *m, double *x, int *n, int *nbeta, double *beta, double *z, double *Dinvfornu, double *logdetDinvfornu, int *family_glmm, double *Dstarinv, double *logdetDstarinv, double *ustar, double *Sigmuhinv, double *logdetSigmuhinv, double *pee, int *nps,  double *b)
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
/*	double *diffs=Calloc(*nps,double);*/
	double *lfuval=Calloc(*m,double);
	double *lfyuval=Calloc(*m,double);
	double *lfutwid=Calloc(*m,double);
/*	double *tempmaxvec=Calloc(*nps,double);*/

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
/*		memset(tempmaxvec,tempmax,*nps);*/
/*		subvec(lfutwidpieces,tempmaxvec,nps,diffs); diffs=lfutwidpieces-tempmaxvec */

		/* calculate diffs= lfutwidpieces-tempmax */
		/* and lfutwid = tempmax + sum(pee[i]*exp(diffs)) */
		for(i=0;i<*nps;i++){
			diffs=lfutwidpieces[i]- tempmax;
			lfutwid[k]+=*(pee+i)*exp(diffs);
		}
	lfutwid[k]=log(lfutwid[k])+tempmax;
	b[k]=*(lfuval+k)+*(lfyuval+k)-*(lfutwid+k);


	}

	Free(Uk);
	Free(xbeta);
	Free(zu);
	Free(eta);
	Free(mzeros);
	Free(lfutwidpieces);

	Free(lfuval);
	Free(lfyuval);
	Free(lfutwid);
}


