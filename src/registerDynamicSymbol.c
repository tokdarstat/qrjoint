#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void INIT(double *, double *, double *, int *, double *, int *, double *, int *, double *, double *, double *, int *);
extern void INIT_noX(double *, double *, int *, double *, double *, int *, double *, double *, double *, int *);
extern void BJQR(double *, double *, double *, int *, double *, int *, double *, int *, double *, double *, double *, double *, int *, int *, double *, int *, double *, double *, double *, int *);
extern void BQDE(double *, double *, int *, double *, double *, int *, double *, double *, double *, double *, int *, int *, double *, int *, double *, double *, double *, int *);
extern void DEV(double *, double *, double *, int *, double *, int *, double *, int *, double *,double *, double *, double *, double *, double *, int *);
extern void DEV_noX(double *, double *, int *, double *, double *, int *, double *, double *, double *, double *, double *, double *, int *);
extern void PRED_noX(double *, double *, double *, int *, double *, double *, double *, int *);
    
static const R_CMethodDef CEntries[] = {
  {"INIT", (DL_FUNC) &INIT, 12}, 
  {"INIT_noX", (DL_FUNC) &INIT_noX, 10},
    {"BJQR", (DL_FUNC) &BJQR, 20},
    {"BQDE", (DL_FUNC) &BQDE, 18},
    {"DEV",  (DL_FUNC) &DEV,  15},
    {"DEV_noX", (DL_FUNC) &DEV_noX,  13},
    {"PRED_noX", (DL_FUNC) &PRED_noX,  8},
    {NULL, NULL, 0}
};

void R_init_sbde(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
