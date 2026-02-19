#ifndef INC_ATOMTYPE_H
#define INC_ATOMTYPE_H
#include "ParameterTypes.h"
//# incl ude "CpptrajStdio.h" // DEBUG
/// Hold parameters for a unique atom type
class AtomType {
    void clearElt() { elt_[0] = ' '; elt_[1] = ' '; elt_[2] = '\0'; }
  public:
    /// Atom hybridization types
    enum HybridizationType { SP = 0, SP2, SP3, UNKNOWN_HYBRIDIZATION };
    /// CONSTRUCTOR
    AtomType() : mass_(0.0), polarizability_(0.0), oidx_(-1), hybrid_(UNKNOWN_HYBRIDIZATION), hasLJ_(false), hasLJ14_(false) { clearElt(); }
    /// CONSTRUCTOR - Mass only
    AtomType(double m) : mass_(m), polarizability_(0.0), oidx_(-1), hybrid_(UNKNOWN_HYBRIDIZATION), hasLJ_(false), hasLJ14_(false) { clearElt(); }
    /// CONSTRUCTOR - Mass, polarizability
    AtomType(double m, double p) : mass_(m), polarizability_(p), oidx_(-1), hybrid_(UNKNOWN_HYBRIDIZATION), hasLJ_(false), hasLJ14_(false) { clearElt(); }
    /// CONSTRUCTOR - Radius, well depth, mass, polarizability
    AtomType(double r, double d, double m, double p) : lj_(r, d), mass_(m), polarizability_(p), oidx_(-1), hybrid_(UNKNOWN_HYBRIDIZATION), hasLJ_(true), hasLJ14_(false) { clearElt(); }
    /// Set type index
    void SetTypeIdx(int i) { oidx_ = i; }
    /// \return default LJ parameters
    LJparmType const& LJ()  const { return lj_; }
    /// \return LJ 1-4 parameters
    LJparmType const& LJ14() const { return lj14_; }
    /// \return True if LJ parameters have been set
    bool HasLJ() const { return hasLJ_; }
    /// \return True if LJ 1-4 parameters have been set
    bool HasLJ14() const { return hasLJ14_; }
    /// \return Atom mass in amu
    double Mass()           const { return mass_;   }
    /// \return Atomic polarizability in Ang^3
    double Polarizability() const { return polarizability_; }
    /// \return Original atom type index. Useful when checking for off-diagonal NB parameters.
    int OriginalIdx()       const { return oidx_; }
    /// \return Atom hybridization
    HybridizationType Hybridization() const { return hybrid_; }
    /// \return Atom type element string
    const char* EltStr() const { return elt_; }
    /// \return true if mass, polarizability, or LJ params are less than incoming
    bool operator<(AtomType const& rhs) const {
      if (FEQ(mass_, rhs.mass_)) {
        if (FEQ(polarizability_, rhs.polarizability_)) {
          if (lj_ == rhs.lj_) {
            return lj14_ < rhs.lj14_;
          } else {
            return lj_ < rhs.lj_;
          }
        } else {
          return polarizability_ < rhs.polarizability_;
        }
      } else {
        return mass_ < rhs.mass_;
      }
    }
    /// \return true if mass, polarizability, and LJ params are the same
    bool operator==(AtomType const& rhs) const {
      return (FEQ(mass_, rhs.mass_) &&
              FEQ(polarizability_, rhs.polarizability_) &&
              lj_ == rhs.lj_ &&
              lj14_ == rhs.lj14_);
    }
    /// Set LJ params
    void SetLJ(LJparmType const& lj) { lj_ = lj; hasLJ_ = true; }
    /// Set LJ 1-4 parameters
    void SetLJ14(LJparmType const& lj) { lj14_ = lj; hasLJ14_ = true; }
    /// Set atom hybridization
    void SetHybridization(HybridizationType h) { hybrid_ = h; }
    /// Set atom type element string
    void SetEltStr(const char* elt) { elt_[0] = elt[0]; elt_[1] = elt[1]; }
    /// \return data size  (2 double for LJparmType)
    static size_t DataSize() { return (4*sizeof(double)) + sizeof(int); }
    /// Assign atom type. Do not overwrite LJ params/hybridization if incoming type does not have them set.
    AtomType& operator=(AtomType const& rhs) {
      if (this == &rhs) return *this;
      //mprintf("DEBUG: Assign elt %c%c oldHasLj=%i oldR=%f newHasLj=%i newR=%f\n", elt_[0], elt_[1], (int)hasLJ_, lj_.Radius(), (int)rhs.hasLJ_, rhs.lj_.Radius());
      if (rhs.hasLJ_)   { lj_   = rhs.lj_;   hasLJ_ = true; }
      if (rhs.hasLJ14_) { lj14_ = rhs.lj14_; hasLJ14_ = true; }
      //mprintf("DEBUG: After Assign elt %c%c newR=%f\n", elt_[0], elt_[1], lj_.Radius());
      mass_ = rhs.mass_;
      polarizability_ = rhs.polarizability_;
      oidx_ = rhs.oidx_;
      if (rhs.hybrid_ != UNKNOWN_HYBRIDIZATION)
        hybrid_ = rhs.hybrid_;
      elt_[0] = rhs.elt_[0];
      elt_[1] = rhs.elt_[1];
      //elt_[2] = rhs.elt_[2];
      //hasLJ_ = rhs.hasLJ_;
      //hasLJ14_ = rhs.hasLJ14_;
      return *this;
    }
    /// COPY CONSTRUCTOR
    AtomType(AtomType const& rhs) :
      lj_(rhs.lj_),
      lj14_(rhs.lj14_),
      mass_(rhs.mass_),
      polarizability_(rhs.polarizability_),
      oidx_(rhs.oidx_),
      hybrid_(rhs.hybrid_),
      hasLJ_(rhs.hasLJ_),
      hasLJ14_(rhs.hasLJ14_)
    {
      elt_[0] = rhs.elt_[0];
      elt_[1] = rhs.elt_[1];
      elt_[2] = '\0';
    }
  private:
    LJparmType lj_;         ///< Default Lennard-Jones parameters (always valid for self).
    LJparmType lj14_;       ///< Lennard-Jones 1-4 parameters (optional).
    double mass_;           ///< Mass in amu
    double polarizability_; ///< Atomic polarizability in Ang^3
    int oidx_;              ///< Original atom type index, for indexing nonbond parameters.
    HybridizationType hybrid_; ///< Atom hybridization
    char elt_[3];              ///< 2 character atom type element string (and null char)
    bool hasLJ_;             ///< True if this type has LJ 1-4 parameters.
    bool hasLJ14_;             ///< True if this type has LJ 1-4 parameters.
};
#endif
