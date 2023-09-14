#ifndef INC_EXEC_CRDTRANSFORM_H
#define INC_EXEC_CRDTRANSFORM_H
#include "Exec.h"
/// Used to transform/condition coordinates 
class Exec_CrdTransform : public Exec {
  public:
    Exec_CrdTransform() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CrdTransform(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int iterativeRmsRefinement(AtomMask const&, bool, double,
                               DataSet_Coords*, DataSet_Coords*) const;
};
#endif
