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

    /// Hold info for each residue
    typedef std::vector<int> Iarray;
    class SymRes {
      public:
        SymRes() {}
        Iarray NonSymAtoms_;
        std::vector<Iarray> SymMasks_;
    };
    std::vector<SymRes> residues_;

    DataSet* rmsd_;
};
#endif
