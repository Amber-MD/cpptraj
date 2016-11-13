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

    bool maponly_;
    Frame* newFrame_;
    Topology* newParm_;
    Topology* stripParm_; // For stripping reference

    Frame rmsRefFrame_;
    Frame rmsTgtFrame_;
    bool rmsfit_;
    DataSet* rmsdata_;
    int tgtPindex_; ///< Topology index of target topology used in mapping
    int tgtNatom_;  ///< Number of atoms in target topology
};
#endif
