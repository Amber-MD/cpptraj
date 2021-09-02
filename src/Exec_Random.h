#ifndef INC_EXEC_RANDOM_H
#define INC_EXEC_RANDOM_H
#include "Exec.h"
/// Set RNG type; create sets with random data. 
class Exec_Random : public Exec {
  public:
    Exec_Random() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Random(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
