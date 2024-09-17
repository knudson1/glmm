#include "myheader.h"
void elc(double *Y, double *X, int *nrowX, int *ncolX, double *eta, int *family,  int *ntrials, double *wts, double *elval, double *elgradient, double *elhessian)
{

	double cumout=0.0;
	double *cpout=R_Calloc(*nrowX, double);
	// unnecessary because calloc zeros memory (malloc doesn't)
	// memset(cpout,0,*nrowX);
	double *cppout=R_Calloc(*nrowX, double);
	// unnecessary because calloc zeros memory (malloc doesn't)
	// memset(cppout,0,*nrowX);

	/*calling cum3, cp3, cpp3 will just change the doubles cumout, cpout, cppout  */
	cum3(eta, nrowX, family, ntrials, wts, &cumout);
	cp3(eta, nrowX, family, ntrials, cpout);
	cpp3(eta, nrowX, family, ntrials, cppout);
    
    /*create diagonal matrix of weights*/
    int sizemat=(*nrowX)*(*nrowX);
    double *wtsmat=R_Calloc(sizemat, double);
    diag(wts,nrowX,wtsmat);
    /*calculate value of Y*wts*/
    int thing1=1,*ione=&thing1;
    double *Ywts=R_Calloc(*nrowX, double);
    matTmatmult(Y,wtsmat,nrowX,ione,nrowX,Ywts);
	/*calculate value of el: Y^T eta-c(eta)  */
	double thing2=0.0, *foo2=&thing2;
	matTvecmult(Ywts,eta,nrowX,ione, foo2); /* Y dot eta goes in foo2 */
	elval[0]=foo2[0]-cumout;
    R_Free(Ywts);
    
    /*calculate Xwts*/
    int size=(*nrowX)*(*ncolX);
    double *Xwts=R_Calloc(size,double);
    matmatmult(wtsmat,X,nrowX,nrowX,ncolX,Xwts);

    R_Free(wtsmat);
	/*calculate gradient of el: X^T (Y-c'(eta))  
	first use loop to calculate Y-c'(eta) and turn c''(eta) into -c''(eta)*/
	double *Yminuscp=R_Calloc(*nrowX,double);
	int i=0;
	for(i=0;i<(*nrowX);i++)
		{
		*(Yminuscp+i)=*(Y+i)-*(cpout+i);
		*(cppout+i)=-*(cppout+i);
		}
	R_Free(cpout);
	/*invoking matTvecmult clobbers dummy values of elgradient with actual values*/
	matTvecmult(Xwts,Yminuscp,nrowX,ncolX,elgradient);
	R_Free(Yminuscp);

	/*  need to create diagonal matrix negcdub=-c''(eta) from vector cppout  */
	double *negcdub=R_Calloc(sizemat,double);
	// unnecessary because calloc zeros memory (malloc doesn't)
	// memset(negcdub,0,sizemat);
	diag(cppout, nrowX, negcdub); 
	R_Free(cppout);

	/*calculate hessian of el: X^T (-c''(eta)) X  
	first start with (-c''(eta))X
 	(nxn)x(nxp)=nxp=nrowX x ncolX*/
	double *mat=R_Calloc((*nrowX)*(*ncolX),double);
	matmatmult(negcdub,X,nrowX,nrowX,ncolX,mat);
	R_Free(negcdub);

	/*then calculate X^T mat    
	this multiplication will clobber elhessian with the correct value*/
	elhessian[0]=1.0;
	matTmatmult(Xwts,mat,nrowX,ncolX,ncolX,elhessian);
	R_Free(mat);
    R_Free(Xwts);
}


