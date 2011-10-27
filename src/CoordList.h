#ifndef INC_COORDLIST_H
#define INC_COORDLIST_H
/// Class: CoordList
/// Hold only the coordinates of frames in an array. Intended for use where
/// complete parm information about frames is not needed and/or memory
/// is an issue. Each frame can have a different number of atoms.
#include "AtomMask.h"
#include <vector>
class CoordList {
    std::vector<float*> coordList;
    std::vector<int> natomList;
    int ncoords;
  public:
    CoordList();
    ~CoordList();

    int AddCoordsByMask(double *, AtomMask *);
    float *Coord(int, int*);

    int Ncoords() { return ncoords; }
};
#endif
