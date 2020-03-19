#ifndef C_WORHP_AUX_H
#define C_WORHP_AUX_H

#include "C_Worhp_Data.h"
#include "C_filter_functions.h"

DLL_PUBLIC void WorhpDiag(OptVar*, Workspace*, Params*, Control*);

DLL_PRIVATE int InitMatrix(Workspace*, WorhpMatrix*, char *s, int s_len);
DLL_PRIVATE int InitMatrixExtend(Workspace*, WorhpMatrix*, char *s, int s_len, int extend);

DLL_PRIVATE void Internal_Init(OptVar*, Workspace*, Params*, Control*, int);

DLL_PUBLIC void PrintIWMT(Workspace*);
DLL_PUBLIC void PrintRWMT(Workspace*);

#endif
