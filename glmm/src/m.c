
#include <R.h>
#include <R_ext/BLAS.h>
#include "myheader.h"

/*  matrix times vector 
a is the matrix (as double)
b is vector (as double)
nrow is number of rows of matrix
ncol is number of cols of matrix
result is where result will be written*/
void matvecmult(double *a, double *b, int *nrow, int *ncol, double *result)
{
    double zero = 0.0;
    double one = 1.0;
    int ione = 1;
    F77_CALL(dgemv)("n", nrow, ncol, &one, a, nrow, b, &ione, &zero, result,
        &ione);
}

/*  transpose(matrix) times vector 
a is the matrix (as double)
b is vector (as double)
nrow is number of rows of matrix BEFORE transpose
ncol is number of cols of matrix BEFORE transpose
result is where result will be written*/
void matTvecmult(double *a, double *b, int *nrow, int *ncol, double *result)
{
    double zero = 0.0;
    double one = 1.0;
    int ione = 1;
    F77_CALL(dgemv)("T", nrow, ncol, &one, a, nrow, b, &ione, &zero, result,
        &ione);
}

/*  matrix times another matrix */
void matmatmult(double *a, double *b, int *nrowa, int *ncola, int *ncolb,
    double *c)
{
    double one = 1.0;
    double zero = 0.0;
    F77_CALL(dgemm)("n", "n", nrowa, ncolb, ncola, &one, a, nrowa, b, ncola,
        &zero, c, nrowa);
}

/*  TRANSPOSED matrix times another matrix 
a is first matrix before transpose
b is second matrix
nrowa number of rows of a BEFORE transpose
ncola number of cols of a BEFORE transpose
ncolb number of cols of b
c to be clobbered to contain answer*/
void matTmatmult(double *a, double *b, int *nrowa, int *ncola, int *ncolb,
    double *c)
{
    double one = 1.0;
    double zero = 0.0;
    F77_CALL(dgemm)("T", "n", ncola, ncolb, nrowa, &one, a, nrowa, b, nrowa,
        &zero, c, ncola);
}

/*  create diagonal matrix from vector */
void diag(double *invec, int *nrow, double *outmat)
{
int i;
int j;
int k=0;
// clang complains about signed to unsigned conversion
unsigned int matsize=(*nrow)*(*nrow);
memset(outmat,0,matsize);

for(i=0;i<(*nrow);++i)
	{
	for(j=0;j<(*nrow);++j)
		{
		if(i==j)		
			outmat[k]=invec[i];
		k++;
		}
	}
}

/*  create identity matrix */
void identmat(int *nrow, int *diagmat)
{
int i;
int j;
int k=0;
// clang complains about signed to unsigned conversion
unsigned int matsize=(*nrow)*(*nrow);
memset(diagmat,0,matsize);

for(i=0;i<(*nrow);++i)
	{
	for(j=0;j<(*nrow);++j)
		{
		if(i==j)		
			diagmat[k]=1;
		k++;
		}
	}
}

/* sum a one-dim array   */
void sumup(double *a, int *lena, double *suma)
{
	*suma=0;	
	int i=0;
	for(i=0;i<(*lena);i++){
		*suma+=*(a+i);
	}
}

/* subtract two vectors: out= a-b */
void subvec(double *a, double *b, int *len, double *out)
{
	int i=0;
	for(i=0;i<(*len);i++){
		*(out+i)=*(a+i)-*(b+i);
	}
}

/* add two vectors: out= a+b */
void addvec(double *a, double *b, int *len, double *out)
{
	int i=0;
	for(i=0;i<(*len);i++){
		*(out+i)=*(a+i)+*(b+i);
	}
}

/* divide one vector by another elementwise: out= a/b */
void divvec(double *a, double *b, int *len, double *out)
{
	int i=0;
	for(i=0;i<(*len);i++){
		*(out+i)=*(a+i)/(*(b+i));
	}
}

/* dot product two vectors or t(a)%*%b */
/*void dot(double *a, double *b, int *len, double *out)*/
/*{*/
/*    double one = 1.0;*/
/*    F77_CALL(ddot)(len,a,one,b,one);*/
/*}*/

/* calculate eigenvalues of a matrix with lower triangle a
the first "n" is to say "just eigvals, not eigvecs"
second argument says lower tri, not upper*/
/*void eigval(double *a, int *ncol, double *result)*/
/*{*/

/*    F77_CALL(dsyev)("N", "L", ncol, a, ncol, b, &ione, &zero, result,*/
/*        &ione);*/
/*}*/

