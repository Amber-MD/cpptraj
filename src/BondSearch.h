#ifndef INC_BONDSEARCH_H
#define INC_BONDSEARCH_H
#include "Topology.h"
/// Search for bonds by distance in given Frame, add to given Topology.
int BondSearch(Topology&, Frame const&, double, int);
int BondSearch_PL(Topology&, Frame const&, double, int);
#endif
