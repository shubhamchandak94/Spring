#ifndef SPRING_QVZ_UTIL_H_
#define SPRING_QVZ_UTIL_H_
/**
 * Utility functions to help do stuff and manage cross-platform issues
 */

#define _CRT_SECURE_NO_WARNINGS

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

//#ifdef LINUX
#include <time.h>
#define _stat stat
#define _alloca alloca
#define restrict __restrict__
//#elif __APPLE__
//#include <time.h>
//#define _stat stat
//#define _alloca alloca
//#else
//#include <malloc.h>
//#include <windows.h>
//#define restrict __restrict
//#endif

namespace spring {
namespace qvz {

// ceiling(log2()) function used in bit calculations
int cb_log2(int x);

// Missing log2 function
#ifndef LINUX
#define log2(x) (log(x) / log(2.0))
#endif

// Missing math symbols
#ifndef INFINITY
#define INFINITY (DBL_MAX + DBL_MAX)
#endif
#ifndef NAN
#define NAN (INFINITY - INFINITY)
#endif
} // namespace qvz
} // namespace spring

#endif
