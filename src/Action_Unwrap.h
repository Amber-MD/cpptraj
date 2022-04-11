#ifndef INC_ACTION_UNWRAP_H
#define INC_ACTION_UNWRAP_H
#include "Action.h"
#include "ImageTypes.h"
namespace Image {
  class List;
}
class Action_Unwrap : public Action {
  public:
    Action_Unwrap();
    ~Action_Unwrap();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Unwrap(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Image::List* imageList_;     ///< List of entities to unwrap
    Image::Mode imageMode_;      ///< Specifiy unwrapping by atom, res, or mol
    std::string maskExpression_; ///< Expression selecting entities to unwrap
    Frame RefFrame_;             ///< Reference frame, updated each DoAction
    Topology* RefParm_;          ///< Reference topology
    bool center_;                ///< If true, determine distances to centers
    Unit allEntities_;           ///< Hold atoms to copy from target to reference
#   ifdef MPI
    Parallel::Comm trajComm_;
#   endif
};
#endif
