#ifndef INC_EXEC_COMPARETOP_H
#define INC_EXEC_COMPARETOP_H
#include "Exec.h"
/// Compare atoms/parameters and report differences between two topologies.
class Exec_CompareTop : public Exec {
  public:
    Exec_CompareTop() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CompareTop(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    void CompareAtoms(Topology const&, Topology const&, CpptrajFile&) const;
    static inline bool Check(bool,bool,const char*, const char*, const char*);
};
#endif
