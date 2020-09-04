#ifndef INC_UNIT_H
#define INC_UNIT_H
#include "Segment.h"
#include <vector>
namespace Cpptraj {
/// Hold 1 or more contiguous segments of atoms
class Unit {
  public:
    Unit() {}

    typedef std::vector<Segment>::const_iterator const_iterator;
    const_iterator segBegin() const { return segArray_.begin(); }
    const_iterator segEnd()   const { return segArray_.end(); }
    unsigned int nSegments()  const {
  private:
    typedef std::vector<Segment> SegArrayType;
    SegArrayType segArray_; ///< Array of segments that make up the Unit
};
}
#endif
