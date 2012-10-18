#ifndef INC_ACTION_PUCKER_H
#define INC_ACTION_PUCKER_H
// Class: Action_Pucker
/// Calculate the ring pucker given 5 atom masks.
#include "Action.h"
class Action_Pucker: public Action {
  public:
    Action_Pucker();
    void print() {}
  private:
    DataSet *puck_;
    AtomMask M1_;
    AtomMask M2_;
    AtomMask M3_;
    AtomMask M4_;
    AtomMask M5_;
    enum PmethodType { ALTONA=0, CREMER };
    PmethodType puckerMethod_;
    bool amplitude_;
    bool useMass_;
    double offset_;
    double puckermin_;
    double puckermax_;

    int init();
    int setup();
    int action();
};
#endif  
