#ifndef INC_EXEC_GRAFT_H
#define INC_EXEC_GRAFT_H
#include "Exec.h"
/// Graft part of one COORDS to another COORDS 
class Exec_Graft : public Exec {
  public:
    Exec_Graft() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Graft(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<int> Iarray;

    static int get_bond_atoms(ArgList&, Iarray&, Iarray&, Topology const&, Topology const&);

    static Topology* modify_top(Topology const&, AtomMask const&, Frame&);

    static DataSet_Coords* get_crd(ArgList&, DataSetList const&, const char*, const char*, Frame&, const char*);

    RetType graft_ic(CpptrajState&, ArgList&) const;
    RetType graft_rms(CpptrajState&, ArgList&) const;
};
#endif
