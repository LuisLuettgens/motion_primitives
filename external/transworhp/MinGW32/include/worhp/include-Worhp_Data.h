<+

  AutoGen5 template

  # Keep all C functions here that are not supposed or
  # not possible to be AutoGen-erated.

+>

/* Defined in Worhp_Data.F90 via C-interop  */
DLL_PUBLIC bool GetUserAction(const Control *cnt, int action);
DLL_PUBLIC void DoneUserAction(Control *cnt, int done);
DLL_PUBLIC void AddUserAction(Control *cnt, int add);
DLL_PUBLIC void SetNextStage(Control*, int);
DLL_PUBLIC int  GetCurrentStage(Control*);
DLL_PUBLIC int  GetPreviousStage(Control*, int);

/* Defined in C_Worhp_Data.c  */
DLL_PUBLIC void WorhpVersion(int *major, int *minor, char (*patch)[PATCH_STRING_LENGTH]);
DLL_PUBLIC int CheckWorhpVersion(int major, int minor, const char *patch);
DLL_PRIVATE size_t MinRWS(OptVar*, Workspace*, Params*, Control*);
DLL_PRIVATE size_t MinIWS(OptVar*, Workspace*, Params*, Control*);

