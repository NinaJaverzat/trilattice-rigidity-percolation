/*
gather_pebble
find_path
reverse_path
written by DaniMuzi https://github.com/DaniMuzi
*/

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
#include "pebble_game.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// PEBBLE GAME



bool gather_pebble (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, int start, bool reverse, std::unordered_set<int> visited)
{
    bool found = 0;
    std::vector<int> path;

    found = find_path(np, pebble_graph, &path, &visited, start);
    if (found && reverse) { reverse_path(np, pebble_graph, path); }


    return found;
}

// Try to gather a pebble at the node start
bool check_rigidity (std::vector<int>* marks, std::vector<int>* np,std::vector<std::vector<int>>* pebble_graph, int start, bool reverse, std::unordered_set<int> visited, std::queue<int>* q)
{
    bool found = 0;

    // Check if there's already a pebble at start
    if ( (*np)[start]>0 )
    {
      found=1;
      (*marks)[start]=-1;
      return found;
    }

    std::vector<int> floppy_path; //path found to pebble
    found = find_path (np, pebble_graph, &floppy_path, &visited, start);

    if (found)  // If a pebble is found, mark the nodes encountered during the pebble rearrangement as floppy (ie nodes on the path leading to the pebble)
    {
       if (reverse) { reverse_path(np, pebble_graph, floppy_path); }
       for (std::vector<int>::iterator node_in_path = floppy_path.begin(); node_in_path!=floppy_path.end(); ++node_in_path)
       {
          if( (*marks)[*node_in_path]==1 ) std::cout << "**WARNING, marking rigid node "<<*node_in_path<<" as floppy**\n";
         (*marks)[*node_in_path]=-1;
       }
    }
    else   // If no pebble is found, mark all visited nodes as rigid and add them to the BFS queue
    {
       for( std::unordered_set<int>::iterator node_in_path=visited.begin(); node_in_path!=visited.end(); ++node_in_path ) 
       {
           if( (*marks)[*node_in_path]==-1 ) std::cout << "**WARNING, marking floppy node "<<*node_in_path<<" as rigid**\n";
           (*marks)[*node_in_path]=1;
           (*q).push(*node_in_path);
       }  
    }

    return found;
}

// Find a path using breath first search.
bool find_path (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, std::vector<int>* path_found, std::unordered_set<int>* visited, int start)
{

    int node, next_node;
    std::vector<int> path;
    std::queue<std::pair<int, std::vector<int>>> q;                             // An element is: a node, the path searched to reach it

    q.push({start, {start}});
    visited->insert(start);

    // Breath first search
    while(!q.empty())
    {
        std::pair<int, std::vector<int>> node_path = q.front();
        q.pop();

        node = node_path.first;
        path = node_path.second;
	
        for (int i=0; i<(int)(*pebble_graph)[node].size(); ++i)                 // Check each neighbor of the current node
        {
            next_node = (*pebble_graph)[node][i];

            if (!visited->contains(next_node))                                  // If the node has not been checked yet
            {
                visited->insert(next_node);                                     // Mark it as checked
                path.push_back(next_node);                                      // The path to next_node is the path to node, plus next node itself

                if ((*np)[next_node] > 0)
                {
                    path_found->swap(path);                                     // next_node has a pebble -> we return the path to it
                    return 1;
                } else
                {
                    q.push({next_node, path});                                  // We will check the neighbors of next_node (if we don't find a pebble first)
                    path.pop_back();
                }
            }
        }
    }

    return 0;

}

// Simply switches each bond i->j in the path into j->i; updates the pebble count
void reverse_path (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, std::vector<int> path)
{
    int node, next_node;

    const int l = (int)path.size();

    node = path[l-1];
    -- (*np)[node];                                                             // Last node in the path loses a pebble

    // Reverse all edges in the pebble graph
    for (int i=l-2; i>=0; --i)
    {
      next_node = path[i];
      (*pebble_graph)[node].push_back(next_node);
      (*pebble_graph)[next_node].erase(std::remove((*pebble_graph)[next_node].begin(), (*pebble_graph)[next_node].end(), node), (*pebble_graph)[next_node].end());
      node = next_node;
    }

    ++(*np)[node];                                                              // First node in the path gains a pebble

}












//
