#ifndef INC_COORDLIST_H
#define INC_COORDLIST_H
#include <vector>
#include "AtomMask.h"
// Class: CoordList
/// Hold only the coordinates of frames in an array. 
/** Intended for use where complete parm information about frames is not 
  * needed and/or memory is an issue. Each frame can have a different 
  * number of atoms.
  */
class CoordList {
    std::vector<float*> coordList;
    std::vector<int> natomList;
    int ncoords;
  public:
    CoordList();
    ~CoordList();

    int AddCoordsByMask(double *, AtomMask *);
    float *Coord(int, int*);
    float *operator[](int);

    int Ncoords() { return ncoords; }
};
#endif
