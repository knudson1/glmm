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
    double *wtsmat=Calloc(sizemat, double);
    diag(wts,nrowX,wtsmat);
    /*calculate value of Y*wts*/
    int thing1=1,*ione=&thing1;
    double *Ywts=Calloc(*nrowX, double);
    matTmatmult(Y,wtsmat,nrowX,ione,nrowX,Ywts);

	/*calculate value of el: Y^T eta-c(eta)  */
	int intone=1;
	double thing2=0.0, *foo2=&thing2;
	matTvecmult(Ywts,eta,nrowX,&intone, foo2); /* Y dot eta goes in foo2 */
	elval[0]=foo2[0]-cumout;
    
    Free(Ywts);
    Free(wtsmat);
}


