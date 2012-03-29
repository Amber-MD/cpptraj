#ifndef INC_COORDLIST_H
#define INC_COORDLIST_H
#include <vector>
#include "AtomMask.h"
#include "Frame.h"
// Class: CoordList
/// Hold only the coordinates of frames in an array. 
/** Intended for use where complete parm information about frames is not 
  * needed and/or memory is an issue. Each frame can have a different 
  * number of atoms.
  */
class CoordList {
  public:
    CoordList();

    int AddFrameByMask(Frame&, AtomMask&);
    std::vector<float> &operator[](int);
    int MaxNatom();
    int Ncoords();

  private:
    std::vector< std::vector<float> > coordList_;
};
#endif
