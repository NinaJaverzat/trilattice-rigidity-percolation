#ifndef R_PERC_H_
#define R_PERC_H_

#include "defs.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// PERCOLATION FUNCTIONS

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



void make_rigid_clusters (std::vector<std::vector<int>>* bond_labels, std::vector<int>* trial_bonds, std::unordered_map<int, int>* RCS, std::unordered_map<int, int>* RCS_dist,
                          std::vector<std::vector<int>>* pebble_graph, std::vector<std::vector<int>>* network,std::vector<int>* np, Scalars* scalars);

// PHASE DIAGRAM OF A SINGLE TRIAL
void single_trial(int MM, std::unordered_map<int, int>* RCS, std::unordered_map<int, int>* RCS_dist, std::vector<int>* bonds, 
            std::vector<std::vector<int>>* network,std::vector<int>* np, std::vector<std::vector<int>>* pebble_graph, Scalars* scalars, OrderParam* ROP, OrderParam* CHI, int k, bool save_conf);



// MOVEMENT ON TRIANGULAR LATTICE

inline int moveright(int site, int L) {
    if (site % L == L - 1) return site - L + 1;
    else return site + 1;
}

inline int moveleft(int site, int L) {
    if (site % L == 0) return site + L - 1;
    else return site - 1;
}

inline int moveup(int site, int L) {
    if (site >= L*(L - 1)) return site - L*(L-1);
    else return site + L;
}

inline int movedown(int site, int L) {
    if (site <= L - 1) return site + L*(L-1);
    else return site - L;
}

// moving to the right
inline int dir0(int site, int L) { return moveright(site, L); }

// moving up (geometrically, up and to the right)
inline int dir1(int site, int L) { return moveup(site, L); }

// moving up and to the left (the moves should commute)
inline int dir2(int site, int L) { return moveleft(moveup(site, L), L); }

// moving to the left
inline int dir3(int site, int L) { return moveleft(site, L); }

// moving to the down (geometrically, down and left)
inline int dir4(int site, int L) { return movedown(site, L); }

// moving down and to the right (the moves should commute)
inline int dir5(int site, int L) { return moveright(movedown(site, L), L); }

inline int choosedir(int site, int d, int L)
{
    switch (d)
    {
      case 0:
            return dir0(site, L);
        case 1:
            return dir1(site, L);
        case 2:
            return dir2(site, L);
        case 3:
            return dir3(site, L);
        case 4:
            return dir4(site, L);
        case 5:
            return dir5(site, L);
    }
    return false;
}

// INDEPENDENT EDGE COVRING

inline void cover_edge(std::vector<std::vector<int>>* pebble_graph, std::vector <int>* np, int from, int to)
{
    --(*np)[from];
    (*pebble_graph)[from].push_back(to);
}

inline int sign(int x){   if (abs(x)<=1) return x;   else return ((x < 0) - (x > 0)); }
inline int x_dist(int u, int v, Scalars* scalars) { return sign( u%(scalars->L) - v%(scalars->L) );}
inline int y_dist(int u, int v, Scalars* scalars) { return sign( u/(scalars->L) - v/(scalars->L) );}


#endif /* R_PERC_H_ */
