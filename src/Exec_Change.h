#ifndef INC_EXEC_CHANGE_H
#define INC_EXEC_CHANGE_H
#include "Exec.h"
class TypeNameHolder;
/// Change things in a topology 
class Exec_Change : public Exec {
  public:
    Exec_Change() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Change(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int ChangeResidueName(Topology&, ArgList&) const;
    int ChangeOresNums(Topology&, ArgList&) const;
    int ChangeIcodes(Topology&, ArgList&) const;
    int ChangeChainID(Topology&, ArgList&) const;
    int ChangeAtomName(Topology&, ArgList&) const;
    static inline int Setup1atomMask(AtomMask&, Topology const&, std::string const&);
    static inline int FindBondTypeIdx(Topology const&, BondArray const&, TypeNameHolder const&);
    int AddBond(Topology&, ArgList&) const;
    int RemoveBonds(CpptrajState&, Topology&, ArgList&) const;
};
#endif
