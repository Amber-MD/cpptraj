#ifndef INC_EXEC_CATCRD_H
#define INC_EXEC_CATCRD_H
#include "Exec.h"
/// Concatenate two or more COORDS data sets 
class Exec_CatCrd : public Exec {
  public:
    Exec_CatCrd() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CatCrd(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
