#ifndef C_STD_H
#define C_STD_H

#if defined (_MSC_VER)
#if defined (__cplusplus)
#define _Bool bool
#else
typedef char _Bool;
typedef char bool;
#endif
#define true 1
#define false 0
#define __bool_true_false_are_defined 1
#define snprintf _snprintf
#else
#include <stdbool.h>
#endif

#include <stdlib.h>

#ifndef NULL
#ifdef __cplusplus
#define NULL 0L
#else
#define NULL ((void *)0L)
#endif
#endif

#define STRING_LENGTH 128

#include "macros.h"
#include "C_worhp_io.h"

typedef int mat_int;    /**< Integer used in WorhpMatrix (int/size_t) */
typedef int rwmt_index;	/**< RWMT entry type for Worhp data structures. */
typedef int iwmt_index;	/**< IWMT entry type for Worhp data structures. */
typedef int counter;    /**< Counter type for Worhp data structures. */

/** Global status flags. Keep in sync with std.f90 */
enum {
  OK               = 0,		/**< @see std::OK */
  notImplemented   = -9000,	/**< @see std::notImplemented */
  notEnoughRWS     = -9001,	/**< @see std::notEnoughRWS */
  notEnoughIWS     = -9002,	/**< @see std::notEnoughIWS */
  missingOptArg    = -9003,	/**< @see std::missingOptArg */
  wrongDimension   = -9004,	/**< @see std::wrongDimension */
  requestRWS       = 9001,	/**< @see std::requestRWS */
  requestIWS       = 9002	/**< @see std::requestIWS */
};

/**
 * Stuff for WORHP's memory management.
 * WORHP uses wMalloc, wCalloc, wRealloc and wFree to manage all
 * user-visible non-temporary memory. Be default WORHP uses the
 * <stdlib.h> functions.
 *
 * typedefs provide prototypes (these should mimic <stdlib.h>, with
 * the possible exception of system-specific __whatever decorations).
 *
 * SetWorhp_X_Function allow setting user-defined memory functions.
 *
 * The actual w_X_ functions are used by WORHP and are user-visible
 * so users can interact with WORHP's memory; doing so may cause
 * Strange Things to happen, so you'll want to be really careful.
 */
typedef void* (*worhpMallocFunction) (size_t size);
typedef void* (*worhpCallocFunction) (size_t num, size_t size);
typedef void* (*worhpReallocFunction) (void* ptr, size_t size);
typedef void  (*worhpFreeFunction) (void* ptr);

DLL_PUBLIC void SetWorhpMallocFunction(worhpMallocFunction f);
DLL_PUBLIC void SetWorhpCallocFunction(worhpCallocFunction f);
DLL_PUBLIC void SetWorhpReallocFunction(worhpReallocFunction f);
DLL_PUBLIC void SetWorhpFreeFunction(worhpFreeFunction f);

DLL_PUBLIC void* wMalloc(size_t size);
DLL_PUBLIC void* wCalloc(size_t num, size_t size);
DLL_PUBLIC void* wRealloc(void* ptr, size_t size);
DLL_PUBLIC void  wFree(void* ptr);

/**
 * Runtime "constant" to be assigned the TRUE value used by the current compiler.
 * This "constant" is used to communicate logical values with Fortran.
 * Initialised with a standard value that \em should be recognised as .TRUE.
 * by any decent compiler.
 * @see C_InitBool
 */
DLL_PUBLIC extern bool FORTRAN_TRUE;


/**
 * Runtime "constant" to be assigned the FALSE value used by the current compiler.
 * This "constant" is used to communicate logical values with Fortran.
 * Initialised with a standard value that \em should be recognised as .FALSE.
 * by any decent compiler.
 * @see C_InitBool
 */
DLL_PUBLIC extern bool FORTRAN_FALSE;


/**
 * Init function for global Fortran boolean values.
 * Needed for compiler-independence, since compilers seem to use different
 * values for TRUE and FALSE, and (bitwise) 1 might not be recognised as TRUE.
 * @param[in] valtrue  \c .TRUE._WORHP_BOOL value.
 * @param[in] valfalse \c .FALSE._WORHP_BOOL value.
 * @see Fortran interface: std::interface::InitBool
 */
DLL_PRIVATE void C_InitBool(const _Bool valtrue, const _Bool valfalse);


/**
 * Implementation provided by Fortran function.
 * @see std::meps
 */
DLL_PRIVATE double meps(void);

#endif
