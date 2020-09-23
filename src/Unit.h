#ifndef INC_UNIT_H
#define INC_UNIT_H
#include "Segment.h"
#include <vector>
//namespace Cpptraj {
/// Hold 1 or more contiguous segments of atoms
class Unit {
    typedef std::vector<Segment> SegArrayType;
  public:
    /// CONSTRUCTOR - Set last index to large negative so first AddIndex() call creates a Segment
    Unit() : lastIdx_(-9999) {}
    /// CONSTRUCTOR - Single segment
    Unit(int beg, int end) : segArray_(SegArrayType(1, Segment(beg, end))), lastIdx_(end-1) {}

    typedef std::vector<Segment>::const_iterator const_iterator;
    const_iterator segBegin() const { return segArray_.begin(); }
    const_iterator segEnd()   const { return segArray_.end(); }
    unsigned int nSegments()  const { return segArray_.size(); }
    Segment const& operator[](int idx) const { return segArray_[idx]; }

    /// \return Total unit size
    unsigned int Size() const {
      unsigned int total = 0;
      for (SegArrayType::const_iterator it = segArray_.begin(); it != segArray_.end(); ++it)
        total += it->Size();
      return total;
    }

    /// \return First index in the unit
    int Front() const { return segArray_.front().Begin(); }
    int Back()  const { return segArray_.back().End(); }

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
    SegArrayType segArray_; ///< Array of segments that make up the Unit
    int lastIdx_;           ///< The last index added to the Unit
};
//}
#endif
