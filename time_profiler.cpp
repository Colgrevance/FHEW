#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
// #include <ctime.h>
#include <cstdlib>

#include "time_profiler.h"

namespace profiler {

	unsigned long long rdtsc() {
		unsigned a, d;

		__asm__ volatile("rdtsc" : "=a" (a), "=d" (d));							// What!? I guess it is calling some sort of assembly thing

		return ((unsigned long long) a) | (((unsigned long long) d) << 32);		// What!?
	}
}