#ifndef INC_EXEC_MUTATE_H
#define INC_EXEC_MUTATE_H
#include "Exec.h"
namespace Cpptraj {
namespace Structure {
class Creator;
}
}
/// <Enter description of Exec_Mutate here>
class Exec_Mutate : public Exec {
  public:
    Exec_Mutate() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Mutate(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    RetType doMutate(CpptrajState&, ArgList&, DataSet_Coords*, DataSet_Coords*, Cpptraj::Structure::Creator const&) const;

    int debug_;
};
#endif
