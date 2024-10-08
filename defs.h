#ifndef DEFS_H_
#define DEFS_H_

// ALL LIBS

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <vector>
#include <stack>
#include <cstdlib> // For RAND_MAX
#include <random>

#include <string.h>
#include <signal.h>
#include <assert.h>
#include <complex.h>

#include <ctype.h>

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

struct Scalars
{

	int n, m;                                                                     // Number of active nodes (n) and bonds (m)
	int indep, red;                                                               // Number of independent / redundant bonds

	int L, N, M, T;                                                               // Lattice size, Number of nodes, Number of Bonds, Number of trials

	long double RCSmax, NRC;                                                      // Size of the largest rigid cluster, number of rigid clusters

        int wrap_state;								      // The wrapping state of the configuration: 0 (does not wrap), +1 (horizontal), -1 (vertical), 2 (both)

	unsigned int seed;
	std::mt19937 gen;

};


struct OrderParam
{

	std::vector<long double> y;                                                   // First moment (average over runs of the variable)
	std::vector<long double> y2;                                                  // Second moment (average over runs of variable*variable)

};


#endif /* DEFS_H_ */
