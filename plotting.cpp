
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
#include "rigidity_percolation.h"

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


bool save_config (std::vector<std::vector<int>> bond_labels, Scalars* scalars, int k)
{
  ostringstream fname;
  fname << "./cfgs/rigid_clusters_L" << scalars->L <<"_step_"<<k<<".txt";

  FILE* file = fopen(fname.str().c_str(), "w");
  if (!file)
  {
      std::cerr << "Failed to open file " << fname.str() << std::endl;
      return 1;
  }
  fprintf(file, "%d %d %d %d %d\n",scalars->L,scalars->L,scalars->L,scalars->L,scalars->L);
  //
  int bond,label,u,v,d;
  for(vector<vector<int>>::iterator lbond = bond_labels.begin();lbond!=bond_labels.end();++lbond)
  {
    bond = (*lbond).at(0);
    label = (*lbond).at(1);
    u = bond/3;
    d = bond%3;
    v = choosedir(u,d,scalars->L);
    fprintf(file, "%f %f %f %f %d\n", u%(scalars->L)+( u/(scalars->L) )/2., ( u/(scalars->L) )*sqrt(3)/2., v%(scalars->L)+( v/(scalars->L) )/2., ( v/(scalars->L) )*sqrt(3)/2., label );
  }

  fflush(file);
  fclose(file);
  return 0;
}

void save_network (std::vector<std::vector<int>>* network, Scalars* scalars, int k)
{
  ostringstream fname;
  fname << "./cfgs/network_L" << scalars->L <<"_step_"<<k<<".txt";

  FILE* file = fopen(fname.str().c_str(), "w");
  fprintf(file, "%d %d\n",scalars->L,scalars->L);
  for(int u = 0;u<scalars->N;++u)
  {
    
    for(vector<int>::iterator v=(*network)[u].begin();v!=(*network)[u].end();++v)
    {
      if(*v<u) {fprintf(file,"%d %d\n", u,*v);}
    }
  }
}




//
