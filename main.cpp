/*
Bond-diluted rigidity percolation.

Activate bonds in random order;
Study connectivity percolation using the NZ algorithm;
If NZ says:
- Same cluster -> Pivoting <-> New rigid cluster made of the new bond only;
- Different cluster ->
-- Same rigid cluster -> Overconstaining <-> Redundant bond
-- Different rigid cluster -> Rigidification <-> Independent bond, play the pebble game to build the new rigid cluster
*/

#include <cstdio>
#include <fstream>
#include <iostream>

#include <stack>
#include <vector>

#include <math.h>
#include <time.h>				// Get clock time, to measure total run time
#include <iomanip>
#include <cstdlib>
#include <algorithm>

#include <bits/stdc++.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <initializer_list>

#include "defs.h"
#include "basic_functions.h"
#include "rigidity_percolation.h"
#include "init.h"
#include "plotting.h"

using namespace std;
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


// MAIN FUNCTION

int main(int argc, char* argv[]) {

    if (argc!=4) {
      std::cerr << "Usage: " << argv[0] << " <LATTICE LINEAR SIZE> <NUMBER OF SAMPLES> <NUMBER OF THE RUN>" << std::endl;
      return 1;
    }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// BASIC VARIABLES

    Scalars scalars;
    OrderParam ROP, CHI, Tpb;                                              // Rigidity order parameter, susceptibility

    scalars.L = pow(2, atoi(argv[1]));                                          // Linear size of the lattice
    scalars.T = atoi(argv[2]);                                                  // Number of trials
    const int num = atoi(argv[3]);                                              // Run number
    basic_init(&scalars, &CHI, &ROP, &Tpb);
    double pc = 0.6602;								// Critical bond filling fraction (Jacobs, Thorpe)
    //double p = 0.66;//atoi(argv[2]);						// Set filling fraction
    //int MM = round(p*(scalars.M));						// Number of bonds to be placed

    double pmin = pc-.05;
    double pmax = pc+.05;              // Range of filling fractions
    int p_steps = 10; 								// Number of points in the phase diagram


// RELEVANT VARIABLES
    //std::vector<int> bonds={0,1,5,14,13,33,16};

    std::vector<int> bonds(scalars.M), np(scalars.N);                           // Order of bonds activation, number of pebbles per node
    std::vector<std::vector<int>> network(scalars.N);                      		// Triangular lattice, filled as bonds get activated
    std::vector<std::vector<int>> pebble_graph(scalars.N);                      // Pebble network: pebble_graph[i] = j <=> A directed bond i->j exists


    std::unordered_map<int, int> RCS, RCS_dist;                                     	// Rigid clusters size; Rigid clusters size distribution

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// OUTPUT FILES

    std::ostringstream ROPfname, CHIfname, LOGfname;
    ROPfname << "./res/ROP_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";               // Rigidity order parameter vs filling proba
    CHIfname << "./res/CHI_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";               // Susceptibility vs filling proba


    FILE* myROPfile = fopen(ROPfname.str().c_str(), "w");
    if (!myROPfile)
    {
        std::cerr << "Failed to open file " << ROPfname.str() << std::endl;
        return 1;
    }
    FILE* myCHIfile = fopen(CHIfname.str().c_str(), "w");
    if (!myCHIfile)
    {
        std::cerr << "Failed to open file " << CHIfname.str() << std::endl;
        return 1;
    }
   


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    clock_t tStart = clock();
    std::cout << "L = "<<scalars.L<<" N = "<<scalars.N<<" Seed =  "<<scalars.seed<<"\n";
    std::cout << "////////////////////////////////\n\n";

    for(int k=0; k<p_steps; ++k)
    {
      double p = pmin + (pmax-pmin)*(double)k/p_steps;
      std::cout << "\np = "<< p<<",\t";
      int MM = round(p*(scalars.M));
      std::cout << MM << " bonds out of " << scalars.M << " will be placed\n\n";

      for(int i=1; i<=scalars.T; ++i)
      {
        init( &RCS, &RCS_dist, &bonds, &network, &np, &pebble_graph, &scalars, &ROP, &CHI, p_steps);
        single_trial(MM, &RCS, &RCS_dist, &bonds, &network, &np, &pebble_graph, &scalars, &ROP, &CHI, k, true);
        std::cout<<"The largest rigid cluster has size "<<scalars.RCSmax<<"\n";
        save_network (&network, &scalars, k);
        
      }
      fprintf(myROPfile, "%f %Lf %Lf\n", p,  ROP.y[k]/scalars.T, ROP.y2[k]/scalars.T);
      fprintf(myCHIfile, "%f %Lf %Lf\n", p,  (CHI.y[k])/scalars.T, (CHI.y2[k])/scalars.T);
      fflush(myROPfile);
      fflush(myCHIfile);
      
    }

    clock_t tEnd = clock();
    std::cout << "\nElapsed time "<< (double)(tEnd - tStart)/CLOCKS_PER_SEC<<"\n";


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


    return 0;

}
