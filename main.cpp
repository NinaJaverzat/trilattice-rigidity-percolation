/*
Bond-diluted rigidity percolation.

Activate bonds in random order;
Identify the rigid clusters
Compute standard percolation quantities [percolation strength, susceptibility, probability to percolate]
Output network information that can be used to visualize rigid clusters
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
    OrderParam ROP, CHI, PP;                                                    // Strength of the rigid cluster (order parameter), susceptibility, probability of percolating

    scalars.L = pow(2, atoi(argv[1]));                                          // Linear size of the lattice
    scalars.T = atoi(argv[2]);                                                  // Number of trials
    const int num = atoi(argv[3]);                                              // Run number
    basic_init(&scalars, &CHI, &ROP, &PP);
    double pc = 0.6602;								// Critical bond filling fraction (Jacobs, Thorpe)
    
    double pmin = pc-.025;
    double pmax = pc+.025;                                                       // Range of filling fractions
    int p_steps = 20; 								// Number of points in the phase diagram


    // RELEVANT VARIABLES
    //std::vector<int> bonds={0,1,5,14,13,33,16};				// To specify a desired list of bonds instead

    std::vector<int> bonds(scalars.M), np(scalars.N);                           // Order of bonds activation, number of pebbles per node
    std::vector<std::vector<int>> network(scalars.N);                      		// Triangular lattice, filled as bonds get activated
    std::vector<std::vector<int>> pebble_graph(scalars.N);                      // Pebble network: pebble_graph[i] = j <=> A directed bond i->j exists


    std::unordered_map<int, int> RCS, RCS_dist;                                     	// Rigid clusters size; Rigid clusters size distribution

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// OUTPUT FILES

    std::ostringstream ROPfname, CHIfname, PPfname, LOGfname;
    ROPfname << "./res/ROP_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";               // Rigidity order parameter vs filling proba
    CHIfname << "./res/CHI_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";               // Susceptibility vs filling proba
    PPfname << "./res/PP_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";                 // Proba to percolate vs filling proba
    LOGfname << "./res/LOG_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";


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
   
    FILE* myPPfile = fopen(PPfname.str().c_str(), "w");
    if (!myPPfile)
    {
        std::cerr << "Failed to open file " << PPfname.str() << std::endl;
        return 1;
    }
    FILE* myLOGfile = fopen(LOGfname.str().c_str(), "w");
    if (!myLOGfile)
    {
        std::cerr << "Failed to open file " << LOGfname.str() << std::endl;
        return 1;
    }

    fprintf(myLOGfile, "# Seed for random numers generation: %u\n", scalars.seed);
    fprintf(myLOGfile, "# L = %d, N = %d, M = %d, T = %d\n", scalars.L, scalars.N, scalars.M, scalars.T);
    fprintf(myLOGfile, "# pmin = %f, pmax = %f, p_steps = %d\n", pmin, pmax, p_steps);
    fflush(myLOGfile);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    cout << "////////////////////////////////////\n";
    cout << "L = "<<scalars.L<<" N = "<<scalars.N << " T = " << scalars.T <<"\tSeed =  "<<scalars.seed<<"\n";
    cout << "////////////////////////////////////\n\n";
    clock_t tStart = clock();
    for(int k=0; k<p_steps; ++k)
    {
      double p = pmin + (pmax-pmin)*(double)k/p_steps;
      std::cout << "p = "<< p<<"\n";
      int MM = round(p*(scalars.M));

      int w[4]={0}; // Pv, Ph, Pb
      for(int i=1; i<=scalars.T; ++i)
      {
        init(&RCS, &RCS_dist, &bonds, &network, &np, &pebble_graph, &scalars, &ROP, &CHI, p_steps);
        single_trial(MM, &RCS, &RCS_dist, &bonds, &network, &np, &pebble_graph, &scalars, &ROP, &CHI, k, false);
        
        // Sanity check
        int npebbles=0;//total number of pebbles in the graph
        for(std::vector<int>::iterator it = np.begin();it!=np.end();++it) npebbles+=*it;
        if(2*scalars.N-scalars.indep!=npebbles) cout << "ERROR, 2N-n_indep = "<<2*scalars.N-scalars.indep<<" while n_pebbles = "<<npebbles<<"\n\n";
        std::cout<<scalars.NRC << " rigid clusters; The largest rigid cluster has size "<<scalars.RCSmax<<" and wrap state "<<scalars.wrap_state<<"\n\n";

        //save_network(&network, &scalars, k); // If you wish to save the network information
        w[scalars.wrap_state+1]++;     
      }

      fprintf(myPPfile, "%f %Lf %Lf %Lf\n", p,  (long double)w[0]/scalars.T, (long double)w[2]/scalars.T, (long double)w[3]/scalars.T);
      fprintf(myROPfile, "%f %Lf %Lf\n", p,  ROP.y[k]/scalars.T, ROP.y2[k]/scalars.T);
      fprintf(myCHIfile, "%f %Lf %Lf\n", p,  (CHI.y[k])/scalars.T, (CHI.y2[k])/scalars.T);
      fflush(myROPfile);
      fflush(myCHIfile);
      fflush(myPPfile);
    }

    
    clock_t tEnd = clock();
    fprintf(myLOGfile, "# Time taken: %.2fs\n", (double)(tEnd - tStart)/CLOCKS_PER_SEC);
    std::cout << "# Time taken: "<< (double)(tEnd - tStart)/CLOCKS_PER_SEC<<"s\n";
    fclose(myPPfile);
    fclose(myROPfile);
    fclose(myCHIfile);
    fflush(myLOGfile);
    fclose(myLOGfile);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


    return 0;

}
