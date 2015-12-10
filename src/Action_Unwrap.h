#ifndef INC_ACTION_UNWRAP_H
#define INC_ACTION_UNWRAP_H
#include "Action.h"
#include "ImageTypes.h"
class Action_Unwrap : public Action {
  public:
    Action_Unwrap();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Unwrap(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Image::PairType imageList_;
    Image::Mode imageMode_;
    std::string maskExpression_;
    Frame RefFrame_;
    Topology* RefParm_;
    bool orthogonal_;
    bool center_;
};
#endif
