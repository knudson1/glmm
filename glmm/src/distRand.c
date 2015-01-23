/* tested in tests/testpiecesBH.R and tests/objfunTest.R */
#include "myheader.h"
void distRandGenC(double *SigmaInv, double *logdet, int *nrow, double *uvec, double *mu, double *distRandGenVal)
{
	double *umu=Calloc(*nrow,double);
	/* calculates umu = uvec-mu */
	subvec(uvec,mu,nrow,umu);

	double *blah=Calloc(*nrow,double);
	/* calculates blah=SigmaInv%*%umu  */
	matvecmult(SigmaInv,umu,nrow,nrow,blah);
	double *blah2=Calloc(1,double);
	int thing1=1,*ione=&thing1;
	/*calculates blah2=t(umu)%*%Sigma.inv %*% umu */
	matTvecmult(umu,blah,nrow,ione,blah2);

	Free(blah);
	Free(umu);
	*distRandGenVal=.5*((-*nrow)*log(2*PI)+ (*logdet)-(*blah2));
	Free(blah2);
}
