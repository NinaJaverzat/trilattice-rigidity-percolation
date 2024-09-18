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
#include "plotting.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////





void make_rigid_clusters (std::vector<std::vector<int>>* bond_labels, std::vector<int>* trial_bonds, std::unordered_map<int, int>* RCS, std::unordered_map<int, int>* RCS_dist,
                          std::vector<std::vector<int>>* pebble_graph, std::vector<std::vector<int>>* network,std::vector<int>* np, Scalars* scalars)
{
    int label = 0;
    queue<int> q;
    int node, r1, r2, d, u, v, dd;
    unordered_set<int>::iterator nextnode;
    int bond, wrap;

    
    while ( !(*trial_bonds).empty() )
    {
        bond = (*trial_bonds).front();
        r1 = bond/3;
        d = bond%3;
        r2 = choosedir(r1, d, scalars->L);
        ++label;
        wrap=0;
        std::vector<int> marks(scalars->N,0);// Rigid (+1) or floppy (-1) node marks
        std::vector<int> dx(scalars->N,0); // x and y displacements along the path, of each node wrt the root of the rigid cluster (r1)
        std::vector<int> dy(scalars->N,0);

        // Gather 3 pebbles at r1 and r2
        while((*np)[r1] < 2 && find_pebble(np, pebble_graph, r1, {}, 1)) {;}
        while((*np)[r2] < 2 && find_pebble(np, pebble_graph, r2, {r1}, 1)) {;}
        if ((*np)[r1]!=2 and (*np)[r2]!=1) cout << "ERROR, could not gather 3 pebbles at "<<r1<<" , "<<r2<<" which should not happen\n";

        // Breadth First Search the mutually rigid nodes
        q.push(r1);
        marks[r1]=1;marks[r2]=1;
        while(!q.empty())
        {
           node = q.front(); q.pop();
              
           // Check the nearest neighbors; for each new unmarked node, perform a pebble search and attempt to free one more pebble keeping the 3 others pinned.
           for(vector<int>::iterator next_node = (*network)[node].begin(); next_node != (*network)[node].end(); ++next_node)
           { 
              // If the next node is not marked yet
              if(marks[*next_node] == 0)
              {
                 // Try to gather a pebble, mark the encountered nodes accordingly and add them to the rigid cluster and the BFS queue
                 check_rigidity(&marks, np, pebble_graph,*next_node, 1, {r1,r2}, &q);
                 marks[r1]=1;marks[r2]=1; // this should be useless
              }
              if(marks[*next_node] == 1) // If next_node has been marked rigid, we compute its displacement to the root node r1 along the path (Machta wrapping criterion)
              {      
                 int dnnn_x = x_dist(*next_node,node,scalars);
                 int dnnn_y = y_dist(*next_node,node,scalars);                 
                 if(dx[*next_node]==0 and dy[*next_node]==0) {dx[*next_node] = dx[node]+dnnn_x;dy[*next_node]=dy[node]+dnnn_y;}
                 else
                 {
                    if ( abs(dx[node]-dx[*next_node]) >1  && abs(dy[node] - dy[*next_node])<=1 ) { if (wrap==-1) {wrap=2;}else if(wrap!=2)wrap=1; }
                    if ( abs(dx[node]-dx[*next_node]) <=1  && abs(dy[node] - dy[*next_node])>1 ) { if (wrap==1) {wrap=2;}else if (wrap!=2)wrap=-1; }
                    if ( abs(dx[node]-dx[*next_node]) >1  && abs(dy[node] - dy[*next_node])>1 ) { wrap=2; }
                 }
              }
           }
           marks[node]=1;
        }
        
        // Label the bonds whose end nodes have been marked rigid, and remove them from the list of bonds
        int bondsize=0;
        vector<int>::iterator b = (*trial_bonds).begin();
        while(b!=(*trial_bonds).end())
        {
           u = *b/3;
           dd = *b%3;
           v = choosedir(u, dd, scalars->L);
           if(marks[u]==1 && marks[v]==1)
           {
              (*bond_labels).push_back({*b,label});
              (*trial_bonds).erase(b);
              b--;
              ++bondsize;
           }
        ++b;
        }
        

        if(bondsize>scalars->RCSmax) {scalars->RCSmax=bondsize;scalars->wrap_state=wrap;}
        (*RCS)[label]=bondsize;
        (*RCS_dist)[bondsize-1]+=1;
    
    } 
    scalars->NRC=label; //number of rigid clusters
}



// PHASE DIAGRAM OF A SINGLE TRIAL: from empty to filled lattice

void single_trial(int MM, std::unordered_map<int, int>* RCS, std::unordered_map<int, int>* RCS_dist, std::vector<int>* bonds, 
            std::vector<std::vector<int>>* network,std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, Scalars* scalars, OrderParam* ROP, OrderParam* CHI, int k, bool save_conf)
{

    int b, d, u, v;
    std::vector<int> trial_bonds;
    std::vector<std::vector<int>> bond_labels;
   

    // Place bonds
    for(int i=scalars->M-1; i>=scalars->M-1 - MM+1; --i)                        // Stop when the desired occupation is reached and go to 2nd step: rigid clusters
    {
        // Select bond
        
        b = (*bonds)[i];
        u = b/3;
        d = b%3;
        v = choosedir(u, d, scalars->L);

        // Activate bond and update n, m
        ++scalars->m;
        if ((*network)[u].size() == 0) ++scalars->n;                            // Node u was not yet active -> number of active nodes ++
        if ((*network)[v].size() == 0) ++scalars->n;                            // Node v was not yet active -> number of active nodes ++
        (*network)[u].push_back(v);
        (*network)[v].push_back(u);
        trial_bonds.push_back(b);
        bonds->pop_back();                                                      // Free some memory 

	// Play the pebble game to determine if the new bond is indep or red
	while((*np)[u] < 2 && find_pebble(np, pebble_graph, u, {}, 1)) {;}
        while((*np)[v] < 2 && find_pebble(np, pebble_graph, v, {u}, 1)) {;}
        if( (*np)[u]==2 && (*np)[v]==2 && (scalars->indep < 2*scalars->N - 3) ) {cover_edge(pebble_graph, np, u, v); ++scalars->indep;}
										// An oriented edge is added in pebble_graph
	      else ++scalars->red;							// No edge added in pebble_graph, but we know a red edge exists since u and v are activated



    }
    
    // Identify rigid clusters
    
    bond_labels.clear();
    make_rigid_clusters(&bond_labels, &trial_bonds, RCS, RCS_dist, pebble_graph, network, np, scalars);

    if (save_conf) save_config (bond_labels, scalars, k);

    // Computing the susceptibility for this trial
    double sum=0;
    double sum2=0;
    // Exclude the largest cluster in the sum
    if(scalars->RCSmax>1)
    {
      for(int s=1;s<scalars->RCSmax;++s){
        
        sum2+=(double)(s*s*(*RCS_dist)[s-1])/scalars->M; sum+=(double)(s*(*RCS_dist)[s-1])/scalars->M;}
    }
    else
    {
      sum2=(double)((*RCS_dist)[0])/scalars->M; sum=(double)((*RCS_dist)[0])/scalars->M;
    }
    

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
