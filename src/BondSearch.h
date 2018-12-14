#ifndef INC_BONDSEARCH_H
#define INC_BONDSEARCH_H
#include "Topology.h"
enum BondSearchType { SEARCH_REGULAR = 0, SEARCH_PAIRLIST, SEARCH_GRID, SEARCH_NONE };
/// Search for bonds in given frame using specified search type.
int BondSearch(Topology&, BondSearchType, Frame const&, double, int);
#endif
