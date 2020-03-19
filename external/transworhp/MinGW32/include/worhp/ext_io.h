#ifndef ext_io_h
#define ext_io_h

#include "macros.h"

typedef void (*ext_io_func) (int mode, const char *s);

DLL_PUBLIC void extprintf_(int mode, const char *message);
DLL_PUBLIC void SetOutputFunction(ext_io_func f);

#endif

