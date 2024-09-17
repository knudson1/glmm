/*
  Y is a matrix of dimension (*nrowX,1), that is, a column vector
  X is a matrix of dimension (*nrowX,*ncolX)
*/
#include "myheader.h"
void elGH(double *Y, double *X, int *nrowX, int *ncolX, double *eta, int *family, int *ntrials, double *wts, double *elgradient, double *elhessian)
{
	double *cpout = R_Calloc(*nrowX,double);
	double *cppout = R_Calloc(*nrowX,double);

	/*calling cp3, cpp3 will just change the doubles cpout, cppout  */
	cp3(eta,nrowX,family,ntrials,cpout);
	cpp3(eta,nrowX,family,ntrials,cppout);

    /*create diagonal matrix of weights*/
    int sizemat=(*nrowX)*(*nrowX);
    double *wtsmat=R_Calloc(sizemat, double);
    diag(wts,nrowX,wtsmat);
    
    /*calculate X*wts*/
    int size=(*nrowX)*(*ncolX);
    double *Xwts=R_Calloc(size,double);
    matmatmult(wtsmat,X,nrowX,nrowX,ncolX,Xwts);
    R_Free(wtsmat);

	/*calculate gradient of el: X^T (Y-c'(eta))  
	first use loop to calculate Y-c'(eta). then turn c''(eta) into -c''(eta)*/
	double *Yminuscp=R_Calloc(*nrowX,double);
	int i=0;
	for(i=0;i<(*nrowX);i++)
		{
		Yminuscp[i]=Y[i]-cpout[i];
		cppout[i]=-cppout[i];
		}
	R_Free(cpout);
	/*invoking matTvecmult clobbers dummy values of elgradient with actual values*/
	matTvecmult(Xwts,Yminuscp,nrowX,ncolX,elgradient);
	R_Free(Yminuscp);

	/*  need to create diagonal matrix negcdub=-c''(eta) from vector cppout  */
	double *negcdub=R_Calloc(sizemat,double);
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
	matTmatmult(Xwts,mat,nrowX,ncolX,ncolX,elhessian);
	R_Free(mat);
    R_Free(Xwts);
}


