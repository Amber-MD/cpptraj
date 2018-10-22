#ifndef INC_EXEC_UPDATEPARAMETERS_H
#define INC_EXEC_UPDATEPARAMETERS_H
#include "Exec.h"
#include "ParameterSet.h"
/// Update parameters in a topology with those from a data set. 
class Exec_UpdateParameters : public Exec {
  public:
    Exec_UpdateParameters() : Exec(HIDDEN) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_UpdateParameters(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    ParameterSet GetParameters(Topology const&) const;
    int UpdateBondParams(ParmHolder<BondParmType>&, ParmHolder<BondParmType> const&) const;
    int UpdateParams(Topology&, ParameterSet const&) const;
};
#endif
