#ifndef INC_EXEC_CHANGE_H
#define INC_EXEC_CHANGE_H
#include "Exec.h"
#include "ParameterHolders.h"
/// Change things in a topology 
class Exec_Change : public Exec {
  public:
    Exec_Change() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Change(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int ChangeResidueName(Topology&, ArgList&) const;
    int ChangeAtomName(Topology&, ArgList&) const;
    static inline int Setup1atomMask(AtomMask&, Topology const&, std::string const&);
    static inline int FindBondTypeIdx(Topology const&, BondArray const&, AtomTypeHolder const&);
    int AddBond(Topology&, ArgList&) const;
};
#endif
