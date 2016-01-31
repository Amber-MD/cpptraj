#ifndef INC_EXEC_PERMUTEDIHEDRALS_H
#define INC_EXEC_PERMUTEDIHEDRALS_H
#include "Exec.h"
// NOTE: Formerly Action_DihedralScan
class Exec_PermuteDihedrals : public Exec {
  public:
    Exec_PermuteDihedrals() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_PermuteDihedrals(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
