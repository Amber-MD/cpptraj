#ifndef INC_ATOMTYPE_H
#define INC_ATOMTYPE_H
#include "ParameterTypes.h"
/// Hold parameters for a unique atom type TODO LJ off diagonal
class AtomType {
  public:
    AtomType() : mass_(0.0), oidx_(-1) {}
    AtomType(double r, double d, int o) : lj_(r, d), mass_(0.0), oidx_(o) {} // TODO deprecate
    /// Mass only
    AtomType(double m) : mass_(m), oidx_(-1) {}
    /// Radius, well depth, mass, original type index
    AtomType(double r, double d, double m, int i) : lj_(r, d), mass_(m), oidx_(i) {}
    /// \return default LJ parameters
    LJparmType const& LJ() const { return lj_; }
    /// \return Atom mass in amu
    double Mass()          const { return mass_;   }
    /// \return Original atom type index. Useful when checking for off-diagonal NB parameters.
    int OriginalIdx()      const { return oidx_; }
    /// \return true if LJ params are less than incoming
    bool operator<(AtomType const& rhs) const { return lj_ < rhs.lj_; }
    /// \return true if LJ params are the same
    bool operator==(AtomType const& rhs) const { return lj_ == rhs.lj_; }
    /// Used to modify LJ params
    LJparmType& SetLJ() { return lj_; }
  private:
    LJparmType lj_; ///< Default Lennard-Jones parameters (always valid for self).
    double mass_;   ///< Mass in amu
    int oidx_; ///< Original atom type index. TODO deprecate
};
#endif
