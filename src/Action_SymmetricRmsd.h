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
    class AtomPair {
      public:
        AtomPair() {}
        AtomPair(Iarray Iin) : AtomIndexes_(Iin), AtomValues_(Iin) {}
        Iarray AtomIndexes_; ///< Indexes into AMap
        Iarray AtomValues_;  ///< Actual atom index at position
    };
    typedef std::vector<AtomPair> PairArray;
    std::vector<PairArray> residues_;

    Iarray AMap_;
    DataSet* rmsd_;
    Frame rmsTgtFrame_;
};
#endif
