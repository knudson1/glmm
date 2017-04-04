
#ifndef GLMM_MYHEADER_H
#define GLMM_MYHEADER_H

#define BERNOULLI 1
#define POISSON 2
#define BINOMIAL 3

#include <R.h>
#include <R_ext/BLAS.h>

void matTvecmult(double *a, double *b, int *nrow, int *ncol, double *result);
void diag(double *invec, int *nrow, double *outmat);
void matTmatmult(double *a, double *b, int *nrowa, int *ncola, int *ncolb,
    double *c);
void matmatmult(double *a, double *b, int *nrowa, int *ncola, int *ncolb,
    double *c);
void matvecmult(double *a, double *b, int *nrow, int *ncol, double *result);
void identmat(int *nrow, int *diagmat);
void sumup(double *a, int *lena, double *suma);
void subvec(double *a, double *b, int *len, double *out);
void addvec(double *a, double *b, int *len, double *out);
void divvec(double *a, double *b, int *len, double *out);

void distRandGenC(double *SigmaInv, double *logdet, int *nrow, double *uvec, double *mu, double *distRandGenVal);
void distRand3C(double *nu, double *mu, int *T, int *nrandom, int *meow, double *Uvec, double *drgradient, double *drhessian);

void cum3(double *etain, int *neta, int *typein, int *ntrials, double *cumout);
void cp3(double *etain, int *neta, int *typein, int *ntrials, double *cpout);
void cpp3(double *etain, int *neta, int *typein, int *ntrials, double *cppout);

void elval(double *Y,  int *nrowX, int *ncolX, double *eta, int *family, int *ntrials, double *elval);
void elGH(double *Y, double *X, int *nrowX, int *ncolX, double *eta, int *family,  int *ntrials, double *elgradient, double *elhessian);

void tdist(double *DstarInv,  int *myq, double *uvec, int *zeta, double *tconst, double *logft);

void objfunc(double *y, double *Umat, int *myq, int *m, double *x, int *n, int *nbeta, double *beta, double *z, double *Dinvfornu, double *logdetDinvfornu, int *family_glmm, double *Dstarinv, double *logdetDstarinv, double *ustar, double *Sigmuhinv, double *logdetSigmuhinv, double *pee, int *nps, int *T, int *nrandom, int *meow, double *nu, int *zeta, double *tconst, double *v, int *ntrials, double *value, double *gradient, double *hessian);

void mcsec(double *gamma, double *thing, double *squaretop, double *numsum, double *y, double *Umat, int *myq, int *m, double *x, int *n, int *nbeta, double *beta, double *z, double *Dinvfornu, double *logdetDinvfornu, int *family_glmm, double *Dstarinv, double *logdetDstarinv, double *ustar, double *Sigmuhinv, double *logdetSigmuhinv, double *pee, int *nps, int *T, int *nrandom, int *meow, double *nu, int *zeta, double *tconst, int *ntrials);


R_CMethodDef cMethods[] = {
	{"matTvecmult" , (DL_FUNC) &matTvecmult, 5, {REALSXP, REALSXP, INTSXP, INTSXP, REALSXP} },
	{"diag" , (DL_FUNC) &diag, 3, {REALSXP,  INTSXP, REALSXP} },
	{"matTmatmult" , (DL_FUNC) &matTmatmult, 6, {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP} },
	{"matmatmult" , (DL_FUNC) &matmatmult, 6, {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP} },
	{"matvecmult" , (DL_FUNC) &matvecmult, 5, {REALSXP, REALSXP, INTSXP, INTSXP, REALSXP} },
	{"identmat" , (DL_FUNC) &identmat, 2, { INTSXP, INTSXP} },
	{"sumup" , (DL_FUNC) &sumup, 3, {REALSXP, INTSXP, REALSXP} },
	{"subvec" , (DL_FUNC) &subvec, 4, {REALSXP, REALSXP, INTSXP, REALSXP} },
	{"addvec" , (DL_FUNC) &addvec, 4, {REALSXP, REALSXP, INTSXP, REALSXP} },
	{"divvec" , (DL_FUNC) &divvec, 4, {REALSXP, REALSXP, INTSXP, REALSXP} },

	{"distRandGenC" , (DL_FUNC) &distRandGenC, 6, {REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP} },
	{"distRand3C" , (DL_FUNC) &distRand3C, 8, {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP} },

	{"cum3" , (DL_FUNC) &cum3, 5, {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP} },
	{"cp3" , (DL_FUNC) &cp3, 5, {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP} },
	{"cpp3" , (DL_FUNC) &cpp3, 5, {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP} },

	{"elval" , (DL_FUNC) &elval, 7, {REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP} },
	{"elGH" , (DL_FUNC) &elGH, 9, {REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP} },

	{"tdist" , (DL_FUNC) &tdist, 6, {REALSXP, INTSXP,  REALSXP, INTSXP, REALSXP, REALSXP} },

	{"objfunc" , (DL_FUNC) &objfunc, 30, {REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP} },

	{"mcsec" , (DL_FUNC) &mcsec, 30, {REALSXP, REALSXP,  REALSXP, REALSXP, REALSXP, REALSXP,  INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP} },


	{NULL, NULL, 0}
}


void R_init_myLib(DLL *info) {R_registerRoutines(info cMethods, NULL, NULL, NULL);}



#include <math.h>

#endif /* GLMM_MYHEADER_H */

