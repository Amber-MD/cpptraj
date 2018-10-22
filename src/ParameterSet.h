#ifndef INC_PARAMETERSET_H
#define INC_PARAMETERSET_H
#include "ParameterTypes.h"
#include "AtomTypeArray.h"
#include "ParameterHolders.h"
/// Hold a set of parameters for atom types, bonds, angles, etc.
class ParameterSet {
  public:
    ParameterSet() : hasLJparams_(false) {}

    AtomTypeArray& AT()                { return atomTypes_; }
    ParmHolder<NonbondType>& NB()      { return nbParm_;    }
    ParmHolder<BondParmType>& BP()     { return bondParm_; }
    ParmHolder<AngleParmType>& AP()    { return angleParm_; }
    ParmHolder<BondParmType>& UB()     { return ubParm_; }
    //ParmHolder<DihedralParmType>& DP() { return dihParm_; }
    ParmHolder<DihedralParmType>& IP() { return impParm_; }
    DihedralParmHolder& DP()           { return dihParm_; }

    void SetHasLJparams(bool b) { hasLJparams_ = b; }
    bool HasLJparams() const { return hasLJparams_; }

    AtomTypeArray const& AT()                const { return atomTypes_; }
    ParmHolder<NonbondType> const& NB()      const { return nbParm_;    }
    ParmHolder<BondParmType> const& BP()     const { return bondParm_; }
    ParmHolder<AngleParmType> const& AP()    const { return angleParm_; }
    ParmHolder<BondParmType> const& UB()     const { return ubParm_; }
    //ParmHolder<DihedralParmType> const& DP() const { return dihParm_; }
    ParmHolder<DihedralParmType> const& IP() const { return impParm_; }
    DihedralParmHolder const& DP()           const { return dihParm_; }

    void Debug(const char*) const;
    void Debug() const { return Debug(""); }
  private:
    AtomTypeArray atomTypes_;
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
