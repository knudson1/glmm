#include <stdlib.h> // for NULL
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

// Only need to declare C functions called from R, that is,
//     called from functions in directory R
//         elc
//         mcsec
//         objfunc
//     called from functions in directory tests
//         cp3
//         cpp3
//         cum3
//         distRand3C
//         distRandGenC
//         elGH
//         elval
//         matTmatmult
//         matTvecmult
//         matvecmult
//         subvec
//         sumup
//         tdist
//
// See Sections 5.4 and 6.15 of Writing R Extensions
// Also see R package fooRegister (in git repo https://github.com/cjgeyer/foo)

#include "myheader.h"

static R_NativePrimitiveArgType elc_types[10] = {REALSXP, REALSXP,
    INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP};

static R_NativePrimitiveArgType mcsec_types[30] = {REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP,
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, REALSXP, INTSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType objfunc_types[30] = {REALSXP, REALSXP,
    INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType cp3_types[5] = {REALSXP, INTSXP, INTSXP,
    INTSXP, REALSXP};

static R_NativePrimitiveArgType distRand3C_types[8] = {REALSXP,
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType distRandGenC_types[6] = {REALSXP,
    REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType elGH_types[9] = {REALSXP, REALSXP,
    INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType elval_types[7] = {REALSXP,  INTSXP,
    INTSXP, REALSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType matTmatmult_types[6] = {REALSXP,
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType matTvecmult_types[5] = {REALSXP,
    REALSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType subvec_types[4] = {REALSXP,
    REALSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType sumup_types[3] = {REALSXP,
    INTSXP, REALSXP};

static R_NativePrimitiveArgType tdist_types[6] = {REALSXP,  INTSXP,
    REALSXP, INTSXP, REALSXP, REALSXP};

// types for function cpp3 below not error, same types as cp3
// types for function cum3 below not error, same types as cp3
// types for function matvecmult below not error, same types
//     as matTvecmult

static R_CMethodDef cMethods[] = {
    {"elc", (DL_FUNC) &elc, 10, elc_types},
    {"mcsec", (DL_FUNC) &mcsec, 30, mcsec_types},
    {"objfunc", (DL_FUNC) &objfunc, 30, objfunc_types},
    {"cp3", (DL_FUNC) &cp3, 5, cp3_types},
    {"cpp3", (DL_FUNC) &cpp3, 5, cp3_types},
    {"cum3", (DL_FUNC) &cum3, 5, cp3_types},
    {"distRand3C", (DL_FUNC) &distRand3C, 8, distRand3C_types},
    {"distRandGenC", (DL_FUNC) &distRandGenC, 6, distRandGenC_types},
    {"elGH", (DL_FUNC) &elGH, 9, elGH_types},
    {"elval", (DL_FUNC) &elval, 7, elval_types},
    {"matTmatmult", (DL_FUNC) &matTmatmult, 6, matTmatmult_types},
    {"matTvecmult", (DL_FUNC) &matTvecmult, 5, matTvecmult_types},
    {"matvecmult", (DL_FUNC) &matvecmult, 5, matTvecmult_types},
    {"subvec", (DL_FUNC) &subvec, 4, subvec_types},
    {"sumup", (DL_FUNC) &sumup, 3, sumup_types},
    {"tdist", (DL_FUNC) &tdist, 6, tdist_types},
    {NULL, NULL, 0, NULL}
};

void attribute_visible R_init_glmm(DllInfo *info)
{
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

