#ifndef C_LEQSOL_DATA_H
#define C_LEQSOL_DATA_H

#include "C_std.h"
#include "C_cs.h"

/**
 *  This struct needs to precisely mirror the
 *  LeqsolWorkspace type in leqsol_data.F90.
 *  To make manual and semi-automatic handling easier
 *  + Keep it sorted
 *    1. By type, descending size (double > [type]* > int > Bool)
 *    2. Inside each type, alphabetically by name.
 *  + Exactly one member declaration per line.
 */
typedef struct LeqsolWorkspaceStruct {
  double EPS;
  double ITREF_TOL;
  double ITSOL_TOL;

  double small;
  double u;

  int    *ITSOL_PCOL;
  int    *ITSOL_PROW;
  double *ITSOL_PVAL;
  double *ITSOL_X0;

  int ITSOL_PCOL_ALLOCATED;
  int ITSOL_PROW_ALLOCATED;
  int ITSOL_PVAL_ALLOCATED;
  int ITSOL_X0_ALLOCATED;

  int DIRECT_METHOD;
  int IPRINT;
  int ITREF_MAXITER;
  int ITSOL_MAXITER;
  int ITSOL_METHOD;
  int ITSOL_PNZ;

  int ordering;
  int print_level;
  int nemin;
  int scaling;

  _Bool ITSOL_PERFORMPRECOND;
  _Bool ITSOL_PSYMMETRIC;
  _Bool SCAL;
  _Bool TRYSIMPLE;

  _Bool solve_blas3;
  _Bool solve_mf;

#ifdef WITH_MA97
  _Bool initsymb;
  _Bool initlu;

  double *valsave;
  double *diag;

  void *akeep;
  void *fkeep;
  void *control;
  void *info;
#endif  /*WITH_MA97*/
} LeqsolWorkspace;

/* Add prototypes of all C functions here, even if they are only Fortran
   wrapper functions, to control their visibility.  */

/**
 * InitLeqsolWorkspace initializes the Leqsol Workspace structure with dummy values
 * @param[inout] leqsolws Set variables in @c leqsolws to dummy values.
 * @author       dlw
 */
DLL_PRIVATE void InitLeqsolWorkspace(LeqsolWorkspace *const leqsolws);
DLL_PRIVATE void ReallocateLeqsolWorkspace(LeqsolWorkspace *const leqsolws, int *const NDIM);
DLL_PRIVATE void FreeLeqsolWorkspace(LeqsolWorkspace *const leqsolws);

#ifdef WITH_MA97
DLL_PRIVATE void AllocateDArray(double **const array, int *const NDIM);
DLL_PRIVATE void FreeDArray(double **array);

#endif

#endif
