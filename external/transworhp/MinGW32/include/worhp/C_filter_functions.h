#ifndef C_FILTER_FUNCTIONS_H
#define C_FILTER_FUNCTIONS_H

#include "C_std.h"
#include "C_Worhp_Data.h"

/* Add prototypes of all C functions here, even if they are only Fortran 
   wrapper functions, to control their visibility.  */

DLL_PRIVATE int InitialiseFilter(const OptVar *opt, Workspace *wsp, const Params *par, const Control *cnt);
DLL_PRIVATE int FreeFilter(const OptVar *opt, Workspace *wsp, const Params *par, const Control *cnt);
DLL_PRIVATE int AddEntryToFilter(const OptVar *opt, Workspace *wsp, const Params *par, const Control *cnt); 
DLL_PRIVATE void PrintFilter(const FilterNode *const first);

#endif
