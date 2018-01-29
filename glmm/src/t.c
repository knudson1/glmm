/* tested in tests/testpiecesBH.R and tests/objfunTest.R 
DstarInv is matrix (D^*)^{-1}
myq is q, length of u
uvec is vector or random effects
zeta is the df of t
tconst is the log of the constants for the t dist
logft is where we put the final calculated value
*/
#include "myheader.h"
void tdist(double *DstarInv,  int *myq, double *uvec, int *zeta, double *tconst, double *logft)
{
	double *blah=Calloc(*myq,double);
	/* calculates blah=DstarInv%*%u  */
	matvecmult(DstarInv,uvec,myq,myq,blah);

	double blah2;
	int thing1=1;
	/*calculates blah2=t(u)%*%Dstarinv %*% u */
	matTvecmult(uvec,blah,myq,&thing1,&blah2);
	Free(blah);


	/* calculates inner=1+blah2/df  */
	double inner=1+blah2/(*zeta);
	double outer=log(inner);
	double thing=  *zeta/2.0+*myq/2.0;
	*logft=*tconst-outer*thing;
}
