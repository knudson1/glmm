/*
  Y is a matrix of dimension (*nrowX,1), that is, a column vector
*/
#include "myheader.h"
#if defined(__GNUC__) || defined(__clang__)
void elval(double *Y,  int *nrowX, int *ncolX __attribute__ ((unused)), double *eta, int *family, int *ntrials, double *wts, double *elval)
#else
void elval(double *Y,  int *nrowX, int *ncolX, double *eta, int *family, int *ntrials, double *wts, double *elval)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
	double cumout = 0.0;

	/*calling cum3, cp3, cpp3 will just change the doubles cumout, cpout, cppout  */
	cum3(eta,nrowX,family,ntrials, wts, &cumout);
    
    /*create diagonal matrix of weights*/
    int sizemat=(*nrowX)*(*nrowX);
    double *wtsmat=R_Calloc(sizemat, double);
    diag(wts,nrowX,wtsmat);
    /*calculate value of Y*wts*/
    int ione = 1;
    double *Ywts=R_Calloc(*nrowX, double);
    matTmatmult(Y,wtsmat,nrowX,&ione,nrowX,Ywts);

	/*calculate value of el: Y^T eta-c(eta)  */
	double foo2;
	matTvecmult(Ywts,eta,nrowX,&ione,&foo2); /* Y dot eta goes in foo2 */
	elval[0]=foo2-cumout;
    
    R_Free(Ywts);
    R_Free(wtsmat);
}


