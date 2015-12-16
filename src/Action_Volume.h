#ifndef INC_ACTION_VOLUME_H
#define INC_ACTION_VOLUME_H
#include "Action.h"
#include "ImagedAction.h"
/// Calculate unit cell volume. 
class Action_Volume: public Action {
  public:
    Action_Volume();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Volume(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    ImagedAction image_;
    DataSet *vol_;
    double sum_;
    double sum2_;
    int nframes_; 
};
#endif
