#ifndef INC_ATOMTYPE_H
#define INC_ATOMTYPE_H
#include "ParameterTypes.h"
/// Hold parameters for a unique atom type TODO LJ off diagonal
class AtomType {
  public:
    AtomType() : radius_(0.0), depth_(0.0), mass_(0.0), oidx_(-1) {}
    AtomType(double r, double d, int o) : radius_(r), depth_(d), mass_(0.0), oidx_(o) {}
    AtomType(double m) : radius_(0.0), depth_(0.0), mass_(m), oidx_(-1) {}
    double Radius() const { return radius_; }
    double Depth()  const { return depth_;  }
    double Mass()   const { return mass_;   }
    void SetRadius(double r) { radius_ = r; }
    void SetDepth(double d)  { depth_ = d; }
    int OriginalIdx() const { return oidx_; }
    bool operator<(AtomType const& rhs)  const {
      return ( (radius_ < rhs.radius_) && (depth_ < rhs.depth_) );
    }
    /// \return true if radius and depth are the same
    bool operator==(AtomType const&) const;
    /// Combine LJ params with this and another type using Lorentz-Berthelot rules
    NonbondType Combine_LB(AtomType const&) const;
    /// \return data size
    static size_t DataSize() { return (3*sizeof(double)) + sizeof(int); }
  private:
    double radius_; ///< VDW radius
    double depth_;  ///< LJ well-depth
    double mass_;   ///< Mass in amu
    int oidx_; ///< Original atom type index.
};
#endif
