#ifndef WORHP_MACROS_H
#define WORHP_MACROS_H

#include "worhp_config.h"
#include "worhp_version.h"

/* 
 * This macro adds functionality similar to C assert to Fortran
 * MAKE_STRING is necessary, because VS treats "\" as escape sequence now (/assume:bscc)
 * and Windows paths contain "\", which may be misinterpreted, or cause warnings.
 */
#ifndef NDEBUG
# define ASSERT(expr) IF (.NOT.(expr)) CALL WORHP_ASSERT("expr", __FILE__, __LINE__)
#else
# define ASSERT(expr)
#endif

#ifdef __GNUC__
# define GCC_VERSION (__GNUC__*10000+__GNUC_MINOR__*100+__GNUC_PATCHLEVEL__)
#else
# define GCC_VERSION -1
#endif

/*
 * This kludge is necessary because the gfortran 4.3 fpp does not
 * behave like standard cpp. This is supposed to be fixed in 4.4.
 */
#if (defined(__INTEL_COMPILER))
# define MAKE_STRING(x) #x
#else
# define MAKE_STRING(x) "x"
#endif

/* Standard way of doing this */
#define Q____(x) #x
#define QUOTE(x) Q____(x)


#ifdef MEM_DEBUG
/* 
 * Make C_realloc calls by Fortran actually suitably call C_DebugRealloc,
 * which is an alias defined in std.F90 to DebugRealloc in C_std.c
 */
# define C_realloc(p,n) C_DebugRealloc(p,n,F2CString(__FILE__),__LINE__)
#endif


/*
 * IVF only defines _WIN32, not WIN32.
 */
#ifdef _WIN32
# ifndef WIN32
#  define WIN32
# endif
#endif

/*
 * Macros for controlling library export and C binding [the latter 
 * being necessary, because !DEC$ STDCALL does not like BIND(C)]
 *
 * USAGE:
 * To decorate a routine, define it in the usual way, and add the
 * macros just behind all USE and IMPLICIT lines, where applicable.
 *
 * TODO:
 * Conditionally use !GNU ATTRIBUTES for GCC >= 4.5.
 *
 * FIXED:
 * ifort's fpp needs the (Linux) -allow nofpp_comments or (Windows)
 * /allow:nofpp_comments switches, or it will eat the !DEC$ line.
 */
#ifndef WIN32

/***** Linux/Unix macros *****/

/* Despite their name, these should apply to both Windows and Linux. 
 * Mingw-w64 cross-compiler says these are not supported. */
#define DLL_PUBLIC  __attribute__ ((visibility("default")))
#define DLL_PRIVATE __attribute__ ((visibility("hidden")))

#define WORHP_CBD    BIND(C)
#define WORHP_CNM(x) BIND(C, name = MAKE_STRING(x))
#define WORHP_API(x) 

#else

#define DLL_PUBLIC
#define DLL_PRIVATE

/***** Windows macros (VS + MinGW) *****/

#ifndef USE_STDCALL

/***** Use default calling convention *****/
#define WORHP_CBD    BIND(C)
#define WORHP_CNM(x) BIND(C, name = MAKE_STRING(x))
#ifdef __GNUC__
        /* Filter this away for GCC, because it does not like $  */
#define WORHP_API(x) 
#else
#define WORHP_API(x) !DEC$ ATTRIBUTES DLLEXPORT :: x
#endif

#else

/***** Use stdcall calling convention *****/
#define WORHP_CBD
#define WORHP_CNM(x)
#ifdef __GNUC__
        /* Filter this away for GCC, because it does not like $  */
#define WORHP_API(x) 
#else
#define WORHP_API(x) !DEC$ ATTRIBUTES DLLEXPORT,STDCALL,REFERENCE,DECORATE::x
#endif

#endif
#endif

/*
 * Macros for secure access to workspace partitions
 * via Workspace Management Tables
 */

/* Get whole slices. Preferred way of access */
#define IWS_SLICE(ENTRY) iws(work%IWMT(ENTRY,3):work%IWMT(ENTRY,4))
#define RWS_SLICE(ENTRY) rws(work%RWMT(ENTRY,3):work%RWMT(ENTRY,4))

/* Get size of the specified slice */
#define IWMT_SIZE(ENTRY) INT(work%IWMT(ENTRY,5), KIND=WORHP_INT)
#define RWMT_SIZE(ENTRY) INT(work%RWMT(ENTRY,5), KIND=WORHP_INT)


/*
 * Macros for custom access with 1-indexing: 1...N
 */

/* Get index in WS */
#define IWS_IDX_1(ENTRY) work%IWMT(ENTRY,1)
#define RWS_IDX_1(ENTRY) work%RWMT(ENTRY,1)

/* Get single element in WS */
#define IWS_ELEM_1(ENTRY,IDX) iws(work%IWMT(ENTRY,1)+IDX)
#define RWS_ELEM_1(ENTRY,IDX) rws(work%RWMT(ENTRY,1)+IDX)

/* Get custom range of WS */
#define IWS_RANGE_1(E,IS,IE) iws(work%IWMT(E,1)+IS:work%IWMT(E,1)+IE)
#define RWS_RANGE_1(E,IS,IE) rws(work%RWMT(E,1)+IS:work%RWMT(E,1)+IE)



/*
 * Macros for custom access with 0-indexing: 0...N-1
 */

/* Get index in WS */
#define IWS_IDX_0(ENTRY) work%IWMT(ENTRY,3)
#define RWS_IDX_0(ENTRY) work%RWMT(ENTRY,3)

/* Get single element in WS */
#define IWS_ELEM_0(ENTRY,IDX) iws(work%IWMT(ENTRY,3)+IDX)
#define RWS_ELEM_0(ENTRY,IDX) rws(work%RWMT(ENTRY,3)+IDX)

/* Get custom range of WS */
#define IWS_RANGE_0(E,IS,IE) iws(work%IWMT(E,3)+IS:work%IWMT(E,3)+IE)
#define RWS_RANGE_0(E,IS,IE) rws(work%RWMT(E,3)+IS:work%RWMT(E,3)+IE)


/*
 * Macros for pointer declaration.
 */
#define INT_POINTER  INTEGER(WORHP_INT),  POINTER, DIMENSION(:)
#define REAL_POINTER REAL(WORHP_DOUBLE),  POINTER, DIMENSION(:)
#define BOOL_POINTER LOGICAL(WORHP_BOOL), POINTER, DIMENSION(:)



/************************************************/
/*                                              */
/*                 NEW MACROS                   */
/* These take more arguments, so the workspace  */
/* instance and iws/rws names can be chosen     */
/* freely (the old macros assume them to be     */
/* "work" and "Xws")                            */
/*                                              */
/************************************************/

/* Get whole slices. Preferred way of access */
#define I_S(wsp,idx,iws) iws(wsp%IWMT(idx,3):wsp%IWMT(idx,4))
#define R_S(wsp,idx,rws) rws(wsp%RWMT(idx,3):wsp%RWMT(idx,4))

/* Get size of the specified slice */
#define I_N(wsp,idx) wsp%IWMT(idx,5)
#define R_N(wsp,idx) wsp%RWMT(idx,5)


/*
 * Macros for custom access with 1-indexing: 1...N
 */

/* Get single element in WS */
#define I_E1(wsp,idx,iws,i) iws(wsp%IWMT(idx,1)+i)
#define R_E1(wsp,idx,rws,i) rws(wsp%RWMT(idx,1)+i)

/* Get custom range of WS */
#define I_R1(wsp,idx,iws,i,j) iws(wsp%IWMT(idx,1)+i:wsp%IWMT(idx,1)+j)
#define R_R1(wsp,idx,rws,i,j) rws(wsp%RWMT(idx,1)+i:wsp%RWMT(idx,1)+j)

/* Get index in WS */
#define I_I(wsp,idx) wsp%IWMT(idx,1)
#define R_I(wsp,idx) wsp%RWMT(idx,1)

/* Version Control */
#if (defined(__GFORTRAN__) || defined(__INTEL_COMPILER))
#define CHECK_WORHP_VERSION IF (F_CheckWorhpVersion(WORHP_MAJOR, WORHP_MINOR, WORHP_PATCH) /= 0) STOP
#else
#define CHECK_WORHP_VERSION if(CheckWorhpVersion(WORHP_MAJOR, WORHP_MINOR, WORHP_PATCH)) {exit(EXIT_FAILURE);}
#endif

#endif
