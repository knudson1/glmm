#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void elc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void matdet(void *, void *, void *);
extern void matinv(void *, void *, void *);
extern void matmatmult(void *, void *, void *, void *, void *, void *);
extern void matsmash(void *, void *, void *, void *);
extern void matsolve(void *, void *, void *, void *);
extern void matvecmult(void *, void *, void *, void *, void *);
extern void mcsec(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void objfunc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void matTvecmult(void *, void *, void *, void *, void *);
extern void diag(void *,void *,void *);
extern void matTmatmult(void *, void *, void *, void *, void *, void *);
extern void identmat(void *, void *);
extern void sumup(void *, void *, void *);
extern void subvec(void *, void *, void *, void *);
extern void addvec(void *, void *, void *, void *);
extern void divvec(void *, void *, void *, void *);

extern void distRandGenC(void *,void *,void *,void *, void *,void *);
extern void distRand3C(void *,void *, void *,void *, void *,void *, void *, void *);

extern void cum3(void *,void *,void *,void *, void *);
extern void cp3(void *,void *,void *,void *, void *);
extern void cpp3(void *,void *,void *,void *, void *);

extern void elval(void *,void *,void *,void *,void *,void *,void *);
extern void elGH(void *,void *,void *,void *,void *,void *,void *,void *,void *);

extern void tdist(void *,void *,void *,void *,void *,void *);


static const R_CMethodDef CEntries[] = {
    {"elc",         (DL_FUNC) &elc,         10},
    {"matdet",      (DL_FUNC) &matdet,       3},
    {"matinv",      (DL_FUNC) &matinv,       3},
    {"matmatmult",  (DL_FUNC) &matmatmult,   6},
    {"matsmash",    (DL_FUNC) &matsmash,     4},
    {"matsolve",    (DL_FUNC) &matsolve,     4},
    {"matvecmult",  (DL_FUNC) &matvecmult,   5},
    {"mcsec",       (DL_FUNC) &mcsec,       30},
    {"objfunc",     (DL_FUNC) &objfunc,     30},
    {"matTvecmult", (DL_FUNC) &matTvecmult,  5},
    {"diag",        (DL_FUNC) &diag,         3},
    {"matTmatmult", (DL_FUNC) &matTmatmult,  6},
    {"identmat",    (DL_FUNC) &identmat,     2},
    {"sumup",       (DL_FUNC) &sumup,        3},
    {"subvec",      (DL_FUNC) &subvec,       4},
    {"addvec",      (DL_FUNC) &addvec,       4},
    {"divvec",      (DL_FUNC) &divvec,       4},
    {"distRandGenC",(DL_FUNC) &distRandGenC, 6},
    {"distRand3C",  (DL_FUNC) &distRand3C,   8},
    {"cum3",        (DL_FUNC) &cum3,         5},
    {"cp3",         (DL_FUNC) &cp3,          5},
    {"cpp3",        (DL_FUNC) &cpp3,         5},
    {"elval",       (DL_FUNC) &elval,        7},
    {"elGH",        (DL_FUNC) &elGH,         9},
    {"tdist",       (DL_FUNC) &tdist,        6},
    {NULL, NULL, 0}
};

void R_init_glmm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

