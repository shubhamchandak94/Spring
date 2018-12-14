#ifndef _UTIL_H_
#define _UTIL_H_
/**
 * Utility functions to help do stuff and manage cross-platform issues
 */

#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef LINUX
	#include <time.h>
	#define _stat stat
	#define _alloca alloca
	#define restrict __restrict__
#elif __APPLE__
    #include <time.h>
    #define _stat stat
    #define _alloca alloca
#else
  #include <malloc.h>
	#include <windows.h>
	#define restrict __restrict
#endif

struct hrtimer_t {
#ifdef LINUX
	struct timespec start;
	struct timespec stop;
	struct timespec res;
#elif __APPLE__
	struct timespec start;
	struct timespec stop;
	struct timespec res;
#else
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
	LARGE_INTEGER freq;
#endif
};

// Cross platform timer interface
void start_timer(struct hrtimer_t *timer);
void stop_timer(struct hrtimer_t *timer);
double get_timer_interval(struct hrtimer_t *timer);

// ceiling(log2()) function used in bit calculations
int cb_log2(int x);

// Missing log2 function
#ifndef LINUX
	#define log2(x) (log(x)/log(2.0))
#endif

// Missing math symbols
#ifndef INFINITY
	#define INFINITY (DBL_MAX + DBL_MAX)
#endif
#ifndef NAN
	#define NAN (INFINITY - INFINITY)
#endif

#endif
