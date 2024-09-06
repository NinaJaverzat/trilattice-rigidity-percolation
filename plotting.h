#ifndef PLOTTING_H_
#define PLOTTING_H_

#include "defs.h"
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


void output_graph (std::vector<std::vector<int>>* pebble_graph, Scalars* scalars);

bool save_config (std::vector<std::vector<int>> bond_labels, Scalars* scalars, int k);

void save_network (std::vector<std::vector<int>>* network, Scalars* scalars, int k);


#endif
