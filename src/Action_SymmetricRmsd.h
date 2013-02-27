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

    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> AtomIndexArray;
    /// Array of groups of potentially symmetric atoms
    AtomIndexArray SymmetricAtomIndices_;

    int debug_;
    Iarray AMap_;       /// AMap_[ref] = tgt
    DataSet* rmsd_;     /// Output DataSet
    Frame remapFrame_;  /// Frame that will contained target re-mapped for symmetry
};
#endif
