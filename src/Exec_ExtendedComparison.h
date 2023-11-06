#ifndef INC_EXEC_EXTENDEDCOMPARISON_H
#define INC_EXEC_EXTENDEDCOMPARISON_H
#include "Exec.h"
/// Calculate extended comparison similarity values for a COORDS set 
class Exec_ExtendedComparison : public Exec {
  public:
    Exec_ExtendedComparison() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ExtendedComparison(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
