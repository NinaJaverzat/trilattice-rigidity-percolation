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
#include "rigidity_percolation.h"
#include "pebble_game.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


void make_rigid_clusters (std::vector<std::vector<int>>* RNlabels, std::vector<std::vector<int>>* RBlabels, std::unordered_map<int, int>* RCS, std::unordered_map<int, int>* RCS_dist,
                          std::vector<std::vector<int>>* pebble_graph, std::vector<std::vector<int>>* network,std::vector<int>* np, Scalars* scalars)
{
    int label = 0;
    queue<int> q;
    int node;
    unordered_set<int>::iterator nextnode;

    for (int r1=0;r1<(scalars->N); ++r1)
    {
        if( (int)(*RBlabels)[r1].size()==0 ) {++label; (*RNlabels)[r1].push_back(label);}		 // Isolated node
        else
        {
         for(vector<int>::iterator r2=(*RBlabels)[r1].begin();r2!=(*RBlabels)[r1].end();++r2)//r1, r2 nodes of the root bond
         {
            // If bond r1-r2 is not labelled yet
            if (*r2<=0)
	          {
               ++label;

               std::vector<int> marks(scalars->N,0);                                         // Rigid (+1) or floppy (-1) node marks
               // Gather 3 pebbles at r1 and r2
               while((*np)[r1] < 2 && gather_pebble(np, pebble_graph, r1, 1, {})) {;}
               while((*np)[-*r2] < 1 && gather_pebble(np, pebble_graph, -*r2, 1, {r1})) {;}
               if ((*np)[r1]!=2 and (*np)[-*r2]!=1) cout << "ERROR, could not gather 3 pebbles at "<<r1<<" , "<<*r2<<" which should not happen\n";

               // Mark r1, r2 nodes as rigid
	             marks[r1]=1;marks[-*r2]=1;

               // BFS the mutually rigid nodes
               q.push(r1);
               while(!q.empty())
               {
                 node = q.front(); q.pop();
                  // Check the nearest neighbors; for each new unmarked node, perform a pebble search and attempt to free one more pebble keeping the 3 others pinned.
                  int k = 0;
                  for(vector<int>::iterator next_node = (*network)[node].begin(); next_node != (*network)[node].end(); ++next_node)
                  {
                     // If the next node is not marked yet
                     if(marks[*next_node] == 0)
                     {
                     // Try to gather a pebble, mark the encountered nodes accordingly and add them to the rigid cluster
                     // If next_node is rigid wrt root bond, add it to the queue
                        check_rigidity(label, &marks, RNlabels, np, pebble_graph, *next_node, 1, {r1,-*r2}, &q);
                        marks[r1]=1;marks[-*r2]=1;

                     }
                     ++k;
                  }
               }

               // Label the bonds
               
               int bondsize=0;
               for(int u=0;u<scalars->N;++u)
               {
                 
                 for(int j = 0;j<(int)(*network)[u].size();++j)
                 {
                   int v = (*network)[u].at(j); //this is always positive
                   
                   if(marks[u]==1 && marks[v]==1)
                   {
                     //if u-v already has a positive label
                     if ((*RBlabels)[u].at(j)>0) cout <<"!!ERROR Trying to label "<<u<<" " <<v<<" from "<<(*RBlabels)[u].at(j)<<" to "<<label<<"\n";
                     else {(*RBlabels)[u].at(j)=label;++bondsize;}
                   }
                 }
               }
               bondsize=bondsize/2;
               if(bondsize>scalars->RCSmax) scalars->RCSmax=bondsize;
               (*RCS)[label]=bondsize;
               (*RCS_dist)[bondsize-1]+=1;
            }
          }
        }
    }
}



// PHASE DIAGRAM OF A SINGLE TRIAL: from empty to filled lattice

void RP_single_trial(int MM, std::vector<std::vector<int>>* RNlabels, std::vector<std::vector<int>>* RBlabels, std::unordered_map<int, int>* RCS, std::unordered_map<int, int>* RCS_dist, std::vector<int>* bonds, std::vector<std::vector<int>>* network,
		std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, Scalars* scalars, OrderParam* ROP, OrderParam* CHI, int k)
{

    int b, d, u, v;

    // Place bonds
    for(int i=scalars->M-1; i>=scalars->M-1 - MM+1; --i)                        // Stop when the desired occupation is reached and go to 2nd step: rigid clusters
    {
        const clock_t t0 = clock();

        // Select bond
        b = (*bonds)[i];
        u = b/3;
        d = b%3;
        v = choosedir(u, d, scalars->L);

        // Activate bond and update n, m
        ++scalars->m;
        if ((*network)[u].size() == 0) ++scalars->n;                            // Node u was not yet active -> number of active nodes ++
        if ((*network)[v].size() == 0) ++scalars->n;                            // Node v was not yet active -> number of active nodes ++
        (*network)[u].push_back(v);(*RBlabels)[u].push_back(-v);
        (*network)[v].push_back(u);(*RBlabels)[v].push_back(-u);                // Initialize RBlabels has a copy of network with negative entries.
        bonds->pop_back();                                                      // Free some memory

	      // Play the pebble game to determine if the new bond is indep or red
	      while((*np)[u] < 2 && gather_pebble(np, pebble_graph, u, 1, {})) {;}
        while((*np)[v] < 2 && gather_pebble(np, pebble_graph, v, 1, {u})) {;}
        if( (*np)[u]==2 && (*np)[v]==2 && (scalars->indep < 2*scalars->N - 3) ) {cover_edge(pebble_graph, np, u, v); ++scalars->indep;}
										// An oriented edge is added in pebble_graph
	      else ++scalars->red;							// No edge added in pebble_graph, but we know a red edge exists since u and v are activated


        const clock_t t1 = clock();
        double dt = (double)(t1 - t0)/CLOCKS_PER_SEC;                           // Time to place one bond


    }

    //cout << (scalars->indep) << " indep bonds, " << (scalars->red) << " red bonds\n";

    // Identify rigid clusters
    make_rigid_clusters(RNlabels, RBlabels, RCS, RCS_dist, pebble_graph, network, np, scalars);
    //cout << (int)(*RCS).size() << " rigid clusters, max. size = "<<scalars->RCSmax<<"\n";

    // Computing the susceptibility for this trial
    double sum=0;
    double sum2=0;
    // Exclude the largest cluster in the sum
    //cout <<(int) (*RCS_dist).size()<< " " << scalars->RCSmax << "\n";
    if(scalars->RCSmax>1)
    {
      for(int s=1;s<scalars->RCSmax;++s){
        //cout << s << " " << (*RCS_dist)[s-1]<<"\n";
        sum2+=(double)(s*s*(*RCS_dist)[s-1])/scalars->M; sum+=(double)(s*(*RCS_dist)[s-1])/scalars->M;}
    }
    else
    {
      sum2=(double)((*RCS_dist)[0])/scalars->M; sum=(double)((*RCS_dist)[0])/scalars->M;
    }
    //cout << sum << " " << sum2 << "\n";

    update_OP(ROP, (double)scalars->RCSmax, (double)(scalars->M), k);
    update_OP(CHI, sum2/sum, 1., k);

    /*cout << "Decomposition in rigid clusters: \n";
    cout << "label size\n";
    for (auto x : *RCS)
    cout << x.first << " " << x.second << endl;*/
    /*cout << "RCS_dist\n";
    for (auto x : *RCS_dist)
    cout << x.first+1 << " " << x.second << endl;*/

}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////




//