#ifndef INC_EXEC_HMASSREPARTITION_H
#define INC_EXEC_HMASSREPARTITION_H
#include "Exec.h"
/// Do hydrogen mass repartitioning on a Topology 
class Exec_HmassRepartition : public Exec {
  public:
    Exec_HmassRepartition() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_HmassRepartition(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    static int repartition(Topology&, double, CharMask const&);
};
#endif
