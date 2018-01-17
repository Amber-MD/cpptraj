#ifndef INC_ATOMTYPE_H
#define INC_ATOMTYPE_H
#include "ParameterTypes.h"
/// Hold LJ params for a unique atom type
class AtomType {
  public:
    AtomType() : radius_(0.0), depth_(0.0), oidx_(-1) {}
    AtomType(double r, double d, int o) : radius_(r), depth_(d), oidx_(o) {}
    double Radius() const { return radius_; }
    double Depth()  const { return depth_;  }
    int OriginalIdx() const { return oidx_; }
    bool operator<(AtomType const& rhs)  const {
      return ( (radius_ < rhs.radius_) && (depth_ < rhs.depth_) );
    }
    /// \return true if radius and depth are the same
    bool operator==(AtomType const&) const;
    /// Combine LJ params with this and another type using Lorentz-Berthelot rules
    NonbondType Combine_LB(AtomType const&) const;
  private:
    NameType name_; ///< Atom type name
    double radius_; ///< VDW radius
    double depth_;  ///< LJ well-depth
    int oidx_; ///< Original atom type index.
};
#endif
