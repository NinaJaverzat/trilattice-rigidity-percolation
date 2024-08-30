#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <math.h>				// Basic math functions
#include <time.h>				// Get clock time, to measure total run time
#include <iomanip>
#include <cstdlib>
#include <algorithm>

#include <bits/stdc++.h>

#include <stdlib.h>     //for using the function sleep
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <initializer_list>

#include "defs.h"
#include "basic_functions.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// INITALIZATION

void basic_init(Scalars* scalars, OrderParam* COP, OrderParam* ROP, OrderParam* Tpb)
{

    scalars->N = scalars->L*scalars->L;
    scalars->M = 3*scalars->N;                                                  // Triangular lattice


    std::random_device rd;
    scalars->seed = rd();                                                       // TO GIVE RANDOM SEED
    scalars->seed = 1413661635;
    // scalars->seed = 508269811;                                               // TO GIVE MANUAL SEED
    scalars->gen = std::mt19937(scalars->seed);

    for (int i=0; i<scalars->M; ++i)
    {
        COP->y.push_back(0);
        COP->y2.push_back(0);
        ROP->y.push_back(0);
        ROP->y2.push_back(0);
        Tpb->y.push_back(0);
        Tpb->y2.push_back(0);
    }

}




/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// GENERIC

// Measurement of a generic order parameter. Here we use the percolation strength

void update_OP(OrderParam* OP, double val, double norm, int ind)
{
    OP->y[ind] += val / norm;
    OP->y2[ind] += (val / norm) * (val / norm);

}


















//
