#ifndef INC_ACTION_PUCKER_H
#define INC_ACTION_PUCKER_H
#include "Action.h"
/// Calculate the ring pucker given 5 or 6 atom masks.
class Action_Pucker: public Action {
  public:
    Action_Pucker();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Pucker(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet *pucker_;
    DataSet* amplitude_;
    DataSet* theta_;
    double puckerMin_;
    double puckerMax_;
    double offset_;
    std::vector<AtomMask> Masks_;
    std::vector<Vec3> AX_;
    enum PmethodType { ALTONA=0, CREMER };
    PmethodType puckerMethod_;
    bool useMass_;
};
#endif  
