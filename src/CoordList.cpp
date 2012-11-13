// CoordList
#include <cstddef>
#include <stdexcept>
#include "CoordList.h"

// CONSTRUCTOR
CoordList::CoordList() { }

// CoordList::AddFrameByMask()
/** Given a frame and a corresponding atom mask,store the selected atoms in a 
  * float array. Currently only for use with Integer masks.
  */
int CoordList::AddFrameByMask(Frame const& frameIn, AtomMask const& MaskIn) {
  if (MaskIn.None()) return 1;
  coordList_.push_back( frameIn.ConvertToFloat(MaskIn) );
  return 0;
}

// CoordList::operator[]()
/** Return a reference to the given coord array. 
  */
std::vector<float> &CoordList::operator[](int crd) {
  if (crd<0 || crd >= (int)coordList_.size()) 
    throw std::range_error("CoordList");
  return coordList_[crd];
}

// CoordList::MaxNatom()
/** Return the size of the largest float array in CoordList.
  */
int CoordList::MaxNatom() {
  unsigned int max = 0;
  for (std::vector< std::vector<float> >::iterator farray = coordList_.begin();
                                                   farray != coordList_.end();
                                                   farray++)
  {
    unsigned int fsize = (*farray).size();
    if ( fsize > max )
      max = fsize;
  }
  // NOTE: Check that max is evenly divisible by 3?
  return (int) (max / 3);
}

// CoordList::Ncoords()
int CoordList::Ncoords() {
  return (int)coordList_.size();
}

