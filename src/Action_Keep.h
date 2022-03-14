#ifndef INC_ACTION_KEEP_H
#define INC_ACTION_KEEP_H
#include "Action.h"
class DataSet_string;
/// Keep only specified parts of the system 
class Action_Keep : public Action {
  public:
    Action_Keep();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Keep(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Action::RetType keepBridge(int, ActionFrame&);

    DataSet_string* bridgeData_; ///< Bridging resdiue ID data set
    int nbridge_;                ///< Number of bridging residues to keep
    std::string bridgeResName_;  ///< Bridging residues name

    AtomMask keepMask_;          ///< Mask of atoms to keep
};
#endif
