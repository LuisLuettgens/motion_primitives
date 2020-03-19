#if defined (_MSC_VER)
#include "fintrf.h"
#if defined(__LP64__) || defined(_M_AMD64) || defined(__amd64)
# define API64(N) INT(N,KIND=8)
#else
# define API64(N) N
#endif
#else
#include "../interfaces/matlab/fintrf_worhp.h"
#endif
