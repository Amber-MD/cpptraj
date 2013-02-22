#ifndef INC_SYMMETRICRMSD_H
#define INC_SYMMETRICRMSD_H
#include "Action.h"
#include "ReferenceAction.h"
#include "RmsAction.h"
/// Action to calculate symmetry-corrected RMSD
class Action_SymmetricRmsd : public Action, ReferenceAction, RmsAction {
  public:
    Action_SymmetricRmsd();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_SymmetricRmsd(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* rmsd_;
};
#endif
