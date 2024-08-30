
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
#include "plotting.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void output_graph (std::vector<std::vector<int>>* pebble_graph, std::vector<int>* np, Scalars* scalars)
{
  for(int u = 0;u<scalars->N;++u)
  {
    //if((int)(*pebble_graph)[u].size()>0)
    //{
      cout << u << " {";
    for(vector<int>::iterator v=(*pebble_graph)[u].begin();v!=(*pebble_graph)[u].end();++v)
    {
      cout << *v << " ";
    }
    cout << "} "<<(*np)[u]<<" pebbles\n\n";
  //}
  }
}



bool save_config (std::vector<std::vector<int>>* network, std::vector<std::vector<int>>* RBlabels, Scalars* scalars, string pathname, double p)
{
  ostringstream fname;
  fname << pathname+"rigid_clusters_L" << scalars->L <<"_p" << p <<".txt";

  FILE* file = fopen(fname.str().c_str(), "w");
  if (!file)
  {
      std::cerr << "Failed to open file " << fname.str() << std::endl;
      return 1;
  }
  fprintf(file, "%d %d %d %d %d\n",scalars->L,scalars->L,scalars->L,scalars->L,scalars->L);
  //
  for(int u=0;u<scalars->N;++u)
  {
    for(int j = 0;j<(int)(*network)[u].size();++j)
    {
      int v = (*network)[u].at(j);
      if(u>v)fprintf(file, "%f %f %f %f %d\n", u%(scalars->L)+( u/(scalars->L) )/2., ( u/(scalars->L) )*sqrt(3)/2., v%(scalars->L)+( v/(scalars->L) )/2., ( v/(scalars->L) )*sqrt(3)/2., (int)(*RBlabels)[u].at(j) );
    }
  }
  fflush(file);
  fclose(file);
  return 0;
}




//
