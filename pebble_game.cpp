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


// Perform a Depth First Search for a pebble to bring to the start node 
bool find_pebble (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, int start, std::unordered_set<int> skips, bool reverse)
{
    std::stack <int> path;
    int node, next_node;
    
    path.push(start);   // At the end of the search, path_found is the (reversed) path to pebble or is empty.

    bool visited[(int)(*np).size()] = {}; // Nodes visited during the search
    visited[start] = 1;
    for (auto &skip : skips) visited[skip] = 1; // Nodes to be skipped are marked as visited

    // Depth first search
    while(!(path).empty())
    {
        node = (path).top();

        for (int i=0; i<(int)(*pebble_graph)[node].size(); ++i)                 // Check each neighbor of the current node
        {
            next_node = (*pebble_graph)[node][i];                                
            if( visited[next_node]==0 )                                         // If the node has not been checked yet
            {                        
                visited[next_node]=1;						// Mark it as checked
                (path).push(next_node);  
                if ((*np)[next_node] > 0)					// If a pebble is found
                {
                    if(reverse) { reverse_path(np, pebble_graph, path); }
                    return 1;
                } 
                break;
            }
        }        
        if (node == (path).top() || (path).size()==0) {(path).pop();}
    }
   return 0;
}

// Reverse edge directions along a path leading to a pebble
void reverse_path (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, std::stack<int> path_found)
{
    int node, next_node;

    node = path_found.top();
    -- (*np)[node];                                                             // First node in the path loses a pebble
    path_found.pop();

    // Reverse all edges in the pebble graph
    while( !path_found.empty() )
    {
      next_node = path_found.top();
      (*pebble_graph)[node].push_back(next_node);
      (*pebble_graph)[next_node].erase(std::remove((*pebble_graph)[next_node].begin(), (*pebble_graph)[next_node].end(), node), (*pebble_graph)[next_node].end());
      node = next_node;
      path_found.pop();
    }

    ++(*np)[node];                                                              // Last node in the path gains a pebble
}

// Same as find_pebble, but we record the path, to be able to mark its nodes later on
bool find_path (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, std::stack<int>* full_path, std::unordered_set<int> skips, int start, std::vector<int>* marks)
{

    int node, next_node;
    
    std::stack<int> pebble_path;
    pebble_path.push(start);          // At the end of the search, path_found is the (reversed) path to pebble or is empty.
    (*full_path).push(start);         // Path is not popped and contains the full path
    
    bool visited[(int)(*np).size()] = {};
    visited[start] = 1;
    for (auto &skip : skips) visited[skip] = 1;

    // Depth first search
    while(!pebble_path.empty())
    {
        node = pebble_path.top();
	
        for (int i=0; i<(int)(*pebble_graph)[node].size(); ++i)                 // Check each neighbor of the current node
        {
            next_node = (*pebble_graph)[node][i];

            if (visited[next_node]==0 && (*marks)[next_node]!=1 )              // If the node has not been checked yet and is not already identified as part of a Laman graph
            {
                visited[next_node]=1;                                     // Mark it as checked
                pebble_path.push(next_node);
                (*full_path).push(next_node);                                      

                if ((*np)[next_node] > 0)
                {
                    full_path->swap(pebble_path);
                    return 1;
                } 
                break;    // To get a well defined path without jumps
            }
        }
        
        if (node == pebble_path.top() || pebble_path.size()==0) {pebble_path.pop(); if(!pebble_path.empty()) (*full_path).push( pebble_path.top() );} // Add the node to which we receede
    }

    return 0;

}

// Try to gather a pebble at the node start and mark the encountered nodes according to the result of the search: floppy (rigid) if the search succeeded (failed)
bool check_rigidity (std::vector<int>* marks,std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, int start, bool reverse, std::unordered_set<int> visited, std::queue<int>* q)
{
    bool found = 0;
    int node;
    // Check if there's already a pebble at start
    if ( (*np)[start]>0 )
    {
      found=1;
      (*marks)[start]=-1;
      return found;
    }

    //std::stack<int> floppy_path; //path found to pebble
    std::stack<int> full_path; // the whole visited path
    found = find_path(np, pebble_graph, &full_path, visited, start, marks);

    // If a pebble is found, mark the nodes encountered during the pebble rearrangement as floppy (ie nodes on the path leading to the pebble)
    if (found)
    {
      if (reverse) { reverse_path(np, pebble_graph, full_path); }

      while (! full_path.empty() )
      {
         node = full_path.top();
         if( (*marks)[node]==1 ) std::cout << "**WARNING, marking rigid node "<<node<<" as floppy**\n";
         (*marks)[node]=-1;
         full_path.pop();      
      }
    }
    else
    {
      // If no pebble is found, mark all visited nodes as rigid and add them to the BFS queue
      while (! full_path.empty() )
      {
         node = full_path.top();
         if( (*marks)[node]==-1 ) std::cout << "**WARNING, marking floppy node "<<node<<" as rigid**\n";
         (*marks)[node]=1;
         (*q).push(node);
         full_path.pop();
         
      }
    }

    return found;
}














//
