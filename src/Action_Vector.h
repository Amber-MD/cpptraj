#ifndef INC_ACTION_VECTOR_H
#define INC_ACTION_VECTOR_H
#include "Action.h"
#include "VectorType.h"
class Action_Vector : public Action {
  public:
    Action_Vector();
  private:
    VectorType *Vector_;

    int init();
};
#endif
