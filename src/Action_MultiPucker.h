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

    Cpptraj::Pucker::PuckerSearch puckerSearch_; ///< Used to search for puckers
    std::vector<DataSet*> data_;                 ///< Output DataSets, 1 per pucker
    Range resRange_;                             ///< Residue range to search
    std::string dsetname_;                       ///< Output data set(s) name
    DataFile* outfile_;                          ///< File to write sets to
    DataSetList* masterDSL_;                     ///< Pointer to master DataSetList
};
#endif
