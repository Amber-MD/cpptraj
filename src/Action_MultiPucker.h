#ifndef INC_ACTION_MULTIPUCKER_H
#define INC_ACTION_MULTIPUCKER_H
#include "Action.h"
#include "Pucker_PuckerSearch.h"
#include "Range.h"
/// Automatically detect and calculate puckers within a residue range. 
class Action_MultiPucker : public Action {
  public:
    Action_MultiPucker();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_MultiPucker(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Cpptraj::Pucker::PuckerSearch puckerSearch_;
    Range resRange_;
    std::string dsetname_;
    DataFile* outfile_;
    DataSetList* masterDSL_;
};
#endif
