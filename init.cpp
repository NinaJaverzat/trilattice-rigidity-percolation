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
#include "init.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// INITALIZATION


void RP_init(std::vector<std::vector<int>>* RNlabels, std::vector<std::vector<int>>* RBlabels, std::unordered_map<int, int>* RCS, std::unordered_map<int, int>* RCS_dist, std::vector<int>* bonds, std::vector<std::vector<int>>* network,
          std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, Scalars* scalars, OrderParam* ROP, OrderParam* CHI, int p_steps)
{

    // Clean datastructures after previous run
    bonds->resize(scalars->M);
    ((*ROP).y).resize(p_steps);
    ((*ROP).y2).resize(p_steps);
    ((*CHI).y).resize(p_steps);
    ((*CHI).y2).resize(p_steps);
    //network->clear();
    RCS->clear();
    RCS_dist->clear();

    // Prepare datastructures for next run

    for (int i=0; i<scalars->N; ++i)                                            // Each node:
    {
        (*network)[i].clear();                                                 // does not have lattice neighbors
        (*RBlabels)[i].clear(); 					       // rigid bond labels
        (*RNlabels)[i].clear();
        (*pebble_graph)[i].clear();                                             // does not have "pebble neighbors"
        (*np)[i] = 2;                                                           // has two pebbles
    }


    // Random order of bond activation
    std::iota(bonds->begin(), bonds->begin()+scalars->M, 0);
    std::shuffle(bonds->begin(), bonds->end(), scalars->gen);

    // Scalar variables
    scalars->n = scalars->m = 0;
    scalars->indep = scalars->red = 0;
    scalars->CSmax = scalars->RCSmax = 1;

}
