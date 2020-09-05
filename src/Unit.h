#ifndef INC_UNIT_H
#define INC_UNIT_H
#include "Segment.h"
#include <vector>
//namespace Cpptraj {
/// Hold 1 or more contiguous segments of atoms
class Unit {
  public:
    /// CONSTRUCTOR - Set last index to large negative so first AddIndex() call creates a Segment
    Unit() : lastIdx_(-9999) {}

    typedef std::vector<Segment>::const_iterator const_iterator;
    const_iterator segBegin() const { return segArray_.begin(); }
    const_iterator segEnd()   const { return segArray_.end(); }
    unsigned int nSegments()  const { return segArray_.size(); }

    /// Add index to the unit
    int AddIndex(int idx) {
      int delta = idx - lastIdx_;
      if (delta < 1)
        return 1;
      else if (delta == 1)
        segArray_.back().UpdateEnd( idx );
      else
        segArray_.push_back( Segment(idx) );
      lastIdx_ = idx;
      return 0;
    }
  private:
    typedef std::vector<Segment> SegArrayType;
    SegArrayType segArray_; ///< Array of segments that make up the Unit
    int lastIdx_;           ///< The last index added to the Unit
};
//}
#endif
