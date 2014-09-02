/* tested in tests/distRandtest.R and tests/objfunTest.R */
#include "myheader.h"
void distRand3C(double *nu, double *mu, int *T, int *nrandom, int *meow, double *Uvec, double *drgradient, double *drhessian)
{
	double *umu=Calloc(*T,double);
	int t;
	int i;
	double *drhessvec=Calloc(*T, double);

	for(t=0;t<*T;t++){
		for(i=*(meow+t);i<*(meow+t+1);i++){
			*(umu+t)+=(*(Uvec+i)-*(mu+t))*(*(Uvec+i)-*(mu+t));
		}
		drgradient[t]=umu[t]/(2*(nu[t])*(nu[t]))-nrandom[t]/(2*nu[t]);
		drhessvec[t]=-umu[t]/((nu[t])*(nu[t])*(nu[t]))+nrandom[t]/(2*(nu[t])*(nu[t]));
	}
	Free(umu);
	diag(drhessvec,T,drhessian);
	Free(drhessvec);
}
