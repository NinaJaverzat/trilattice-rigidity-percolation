#ifndef INIT_H_
#define INIT_H_

#include "defs.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// INITALIZATION

void RP_init(std::vector<std::vector<int>>* RNlabels, std::vector<std::vector<int>>* RBlabels, std::unordered_map<int, int>* RCS, std::unordered_map<int, int>* RCS_dist, std::vector<int>* bonds, std::vector<std::vector<int>>* network,
          std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, Scalars* scalars, OrderParam* ROP, OrderParam* CHI, int p_steps);


#endif /* INIT_H_ */
