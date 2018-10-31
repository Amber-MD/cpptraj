#ifndef INC_PARAMETERSET_H
#define INC_PARAMETERSET_H
#include "ParameterTypes.h"
#include "ParameterHolders.h"
#include "AtomType.h"
/// Hold a set of parameters for atom types, bonds, angles, etc.
class ParameterSet {
  public:
    ParameterSet() : hasLJparams_(false) {}

    ParmHolder<AtomType>& AT()         { return atomTypes_; }
    ParmHolder<NonbondType>& NB()      { return nbParm_;    }
    ParmHolder<BondParmType>& BP()     { return bondParm_; }
    ParmHolder<AngleParmType>& AP()    { return angleParm_; }
    ParmHolder<BondParmType>& UB()     { return ubParm_; }
    //ParmHolder<DihedralParmType>& DP() { return dihParm_; }
    ParmHolder<DihedralParmType>& IP() { return impParm_; }
    DihedralParmHolder& DP()           { return dihParm_; }

    void SetHasLJparams(bool b) { hasLJparams_ = b; }
    bool HasLJparams() const { return hasLJparams_; }

    ParmHolder<AtomType> const& AT()         const { return atomTypes_; }
    ParmHolder<NonbondType> const& NB()      const { return nbParm_;    }
    ParmHolder<BondParmType> const& BP()     const { return bondParm_; }
    ParmHolder<AngleParmType> const& AP()    const { return angleParm_; }
    ParmHolder<BondParmType> const& UB()     const { return ubParm_; }
    //ParmHolder<DihedralParmType> const& DP() const { return dihParm_; }
    ParmHolder<DihedralParmType> const& IP() const { return impParm_; }
    DihedralParmHolder const& DP()           const { return dihParm_; }

    void Debug(const char*) const;
    void Debug() const { return Debug(""); }

    /// Used to track what parameters were updated during UpdateParams
    class UpdateCount {
      public:
        UpdateCount() : nBondsUpdated_(0), nAnglesUpdated_(0),
                        nDihedralsUpdated_(0), nImpropersUpdated_(0),
                        nUreyBradleyUpdated_(0), nAtomTypeUpdated_(0),
                        nLJparamsUpdated_(0) {}
        unsigned int nBondsUpdated_;
        unsigned int nAnglesUpdated_;
        unsigned int nDihedralsUpdated_;
        unsigned int nImpropersUpdated_;
        unsigned int nUreyBradleyUpdated_;
        unsigned int nAtomTypeUpdated_;
        unsigned int nLJparamsUpdated_;
    };
    /// Update this set with parameters from given set
    int UpdateParams(ParameterSet const&, UpdateCount&);
  private:
    //AtomTypeArray atomTypes_;
    ParmHolder<AtomType> atomTypes_;
    ParmHolder<NonbondType> nbParm_;
    ParmHolder<BondParmType> bondParm_;
    ParmHolder<AngleParmType> angleParm_;
    ParmHolder<BondParmType> ubParm_;
    //ParmHolder<DihedralParmType> dihParm_;
    ParmHolder<DihedralParmType> impParm_;
    DihedralParmHolder dihParm_;
    bool hasLJparams_;
};
#endif
