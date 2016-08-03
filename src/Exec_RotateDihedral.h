#ifndef INC_EXEC_ROTATEDIHEDRAL_H
#define INC_EXEC_ROTATEDIHEDRAL_H
#include "Exec.h"
/// Rotate a single dihedral
class Exec_RotateDihedral : public Exec {
  public:
    Exec_RotateDihedral() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_RotateDihedral(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
