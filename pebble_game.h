#ifndef PG_H_
#define PG_H_

#include "defs.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// PEBBLE GAME FUNCTIONS

bool find_pebble (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, int start, std::unordered_set<int> skips, bool reverse);

void reverse_path (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, std::stack<int> path_found);

bool find_path (std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, std::stack<int>* full_path, std::unordered_set<int> skips, int start, std::vector<int>* marks);

bool check_rigidity (std::vector<int>* marks,std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, int start, bool reverse, std::unordered_set<int> visited, std::queue<int>* q);


#endif /* PG_H_ */
