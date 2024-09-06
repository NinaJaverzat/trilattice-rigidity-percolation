#ifndef PG_H_
#define PG_H_

#include "defs.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// PEBBLE GAME FUNCTIONS

bool gather_pebble (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, int start, bool reverse, std::unordered_set<int> visited);

bool check_rigidity (std::vector<int>* marks, std::vector<int>* np,std::vector<std::vector<int>>* pebble_graph, int start, bool reverse, std::unordered_set<int> visited, std::queue<int>* q);

bool find_path (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, std::vector<int>* path_found, std::unordered_set<int>* visited, int start);

void reverse_path (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, std::vector<int> path);


#endif /* PG_H_ */
