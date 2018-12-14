#ifndef INC_BONDSEARCH_H
#define INC_BONDSEARCH_H
#include "Topology.h"
enum BondSearchType { REGULAR = 0, PAIRLIST, GRID };
/// Search for bonds by distance in given Frame, add to given Topology.
int BondSearch(Topology&, Frame const&, double, int);
/// Search for bonds by distance in given Frame using a pair list
int BondSearch_PL(Topology&, Frame const&, double, int);
/// Search for bonds in given frame.
int BondSearch(Topology&, BondSearchType, Frame const&, double, int);
#endif
