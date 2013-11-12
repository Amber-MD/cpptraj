#ifndef INC_CRDTYPE_H
#define INC_CRDTYPE_H
#include <vector>
/// Single-precision representation of coordinates.
/** This class is used by DataSet_Coords to hold coordinates in 
  * single-precision for reduced memory footprint. Layout is:
  * CX0, CY0, CZ0, CX1, ... , CZN
  * VX0, VY0, VZ0, VX1, ... , VZN (if hasVelocities_)
  * BX BY BZ BA BB BG             (if numBoxCrd_ > 0)
  */
class CRDtype {
  public:
    CRDtype() : numCoords_(0), numBoxCrd_(0), hasVelocities_(false) {}
    inline CRDtype(const CRDtype&);
    inline CRDtype& operator=(const CRDtype&);
    void reserve(size_t x, size_t v, size_t b)    {
      numCoords_ = x;
      numBoxCrd_ = b;
      hasVelocities_ = (v > 0);
      farray_.reserve( x + v + b );
    }
    void push_back( double d )                    { farray_.push_back( (float)d ); }
    inline float const& operator[](int idx) const { return farray_[idx];           }
    typedef std::vector<float>::const_iterator const_iterator;
    size_t NumCoords()                      const { return numCoords_;      }
    bool HasVelocities()                    const { return hasVelocities_;  }
    const_iterator CrdBegin()               const { return farray_.begin(); }
    const_iterator CrdEnd()                 const { return farray_.begin() + numCoords_; }
    inline const_iterator VelBegin()               const;
    inline const_iterator VelEnd()                 const;
    const_iterator BoxBegin()               const { return farray_.end() - numBoxCrd_; }
    const_iterator BoxEnd()                 const { return farray_.end();              }
  private:
    std::vector<float> farray_; ///< Hold coords and opt. velo./box info.
    size_t numCoords_;          ///< Number of coordinates.
    size_t numBoxCrd_;          ///< Number of box coordinates.
    bool hasVelocities_;        ///< If true, velocities are stored as well.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// COPY CONSTRUCTOR
CRDtype::CRDtype( const CRDtype& rhs) :
  farray_(rhs.farray_), numCoords_(rhs.numCoords_), numBoxCrd_(rhs.numBoxCrd_),
  hasVelocities_(rhs.hasVelocities_) {}
// ASSIGNMENT
CRDtype& CRDtype::operator=(const CRDtype& rhs) {
  if (this == &rhs) return *this;
  farray_ = rhs.farray_;
  numCoords_ = rhs.numCoords_;
  numBoxCrd_ = rhs.numBoxCrd_;
  hasVelocities_ = rhs.hasVelocities_;
  return *this;
}
// CRDtype::VelBegin() 
CRDtype::const_iterator CRDtype::VelBegin() const {
  if (hasVelocities_)
    return CrdEnd();
  else
    return farray_.end();
}
// CRDtype::VelEnd()
CRDtype::const_iterator CRDtype::VelEnd() const {
  if (hasVelocities_)
    return CrdEnd() + numCoords_;
  else
    return farray_.end();
}
#endif
