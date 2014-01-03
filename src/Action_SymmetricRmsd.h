#ifndef INC_SYMMETRICRMSD_H
#define INC_SYMMETRICRMSD_H
#include "Action.h"
#include "ReferenceAction.h"
#include "RmsAction.h"
/// Action to calculate symmetry-corrected RMSD
class Action_SymmetricRmsd : public Action {
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
    bool remap_;          ///< If true remap symmetric atoms
    Iarray AMap_;         ///< AMap_[ref] = tgt
    DataSet* rmsd_;       ///< Output DataSet
    Frame remapFrame_;    ///< Target frame re-mapped for symmetry
    ReferenceAction REF_; ///< Hold reference frame/traj/options
    RmsAction RMS_;       ///< RMSD-related options/actions
};
#endif
