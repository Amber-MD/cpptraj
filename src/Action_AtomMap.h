#ifndef INC_ACTION_ATOMMAP_H
#define INC_ACTION_ATOMMAP_H
#include "Action.h"
/// Action used to map and reorder atoms in target to reference.
class Action_AtomMap : public Action {
  public:
    Action_AtomMap(); 
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_AtomMap(); }
    void Help() const;
    ~Action_AtomMap();
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet_Coords_REF* TgtFrame_;
    DataSet_Coords_REF* RefFrame_;
    int debug_;
    std::vector<int> AMap_;

    bool maponly_;      ///< If true only generate map
    bool byResidue_;    ///< If true create map on residue by residue basis
    Frame* newFrame_;   ///< Frame for re-mapped target
    Topology* newParm_; ///< Topology for re-mapped target

    Frame rmsRefFrame_; ///< Ref frame for calculating RMS of remapped target to ref
    Frame rmsTgtFrame_; ///< Tgt frame for calculating RMS of remapped target to ref
    bool rmsfit_;       ///< If true, attempt to RMS-fit remapped target to ref
    DataSet* rmsdata_;  ///< RMS of remapped target to ref.
};
#endif
