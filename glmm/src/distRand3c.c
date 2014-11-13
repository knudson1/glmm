/* tested in tests/distRandtest.R and tests/objfunTest.R */
#include "myheader.h"
void distRand3C(double *nu, double *mu, int *T, int *nrandom, int *meow, double *Uvec, double *drgradient, double *drhessian)
{
	double *umu=Calloc(*T,double);
	double *drhessvec=Calloc(*T, double);

	for(int t=0;t<*T;t++){
		for(int i=meow[t];i<meow[t+1];i++){
			umu[t]+=(Uvec[i]-mu[t])*(Uvec[i]-mu[t]);
		}
		drgradient[t]=umu[t]/(2*(nu[t])*(nu[t]))-nrandom[t]/(2*nu[t]);
		drhessvec[t]=-umu[t]/((nu[t])*(nu[t])*(nu[t]))+nrandom[t]/(2*(nu[t])*(nu[t]));
	}
	Free(umu);
	diag(drhessvec,T,drhessian);
	Free(drhessvec);
}
