#ifndef INC_ACTION_PUCKER_H
#define INC_ACTION_PUCKER_H
#include "Action.h"
/// Calculate the ring pucker given 5 or 6 atom masks.
class Action_Pucker: public Action {
  public:
    Action_Pucker();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Pucker(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet *pucker_;
    DataSet* amplitude_;
    DataSet* theta_;
    double puckerMin_;
    double puckerMax_;
    double offset_;
    std::vector<AtomMask> Masks_;
    std::vector<Vec3> AX_;
    enum PmethodType { UNSPECIFIED = 0, ALTONA, CREMER };
    PmethodType puckerMethod_;
    bool useMass_;
};
#endif  
