#ifndef INC_ACTION_SURF_H
#define INC_ACTION_SURF_H
// Surf
#include "Action.h"
#include <vector>

class Surf: public Action {
    DataSet *surf;
    AtomMask Mask1;
    double *distances;
    int soluteAtoms;
    std::vector< std::vector<int> > NeighborList;

    double CalcLCPO(int);
  public:
    Surf();
    ~Surf();

    int init();
    int setup();
    int action();
};
#endif
