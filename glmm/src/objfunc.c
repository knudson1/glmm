#include "myheader.h"
/*x is n by nbeta matrix
beta has length nbeta
Umat is myq by m matrix, one COLUMN of Umat used at a time. Umat must be the R mat transposed 
z is n by myq matrix
*/
void objfunc(double *Umat, int *myq, int *m, double *x, int *n, int *nbeta, double *beta, double *z, double *Dinvfornu, double *logdetDinvfornu, double *lfuval, double *out)
{
	double *Uk=Calloc(myq,double);
	int k=0,i=0,Uindex=0;

	double db1=0.0,*double1=&db1; /*temp to hold info that will be put into lfuval*/

	/* calculate xbeta, used to calculate eta for each Uk=U[k,] in R notation */
	double *xbeta=Calloc(n,double);
	matvecmult(x,beta,n,nbeta,xbeta);
	double *zu=Calloc(n,double);
	double *eta=Calloc(n,double);
	double *mzeros=Calloc(m,double);
	memset(mzeros,0,*m);

	for(k=0;k<*m;k++){
		/*start by getting Uk  */
		for(i=0;i<*myq;i++){
			out[i]=Uk[i]=Umat[Uindex];
			Uindex++;
		}

		/* calculate eta for this value of Uk
		first calculate ZUk*/
		matvecmult(z,Uk,n,myq,zu);

		/* then add xbeta+zu */
		addvec(xbeta,zu,n,eta);

		/*log f_theta(u_k)*/
		distRandGenC(Dinvfornu,logdetDinvfornu,myq,Uk,mzeros,double1); 
		*(lfuval+k)=*double1; 
		
	}

	Free(Uk);
	Free(xbeta);
	Free(zu);
	Free(eta);
	Free(mzeros);
}


