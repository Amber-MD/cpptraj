#ifndef INC_ACTION_VECTOR_H
#define INC_ACTION_VECTOR_H
#include "Action.h"
#include "VectorType.h"
class Action_Vector : public Action {
  public:
    Action_Vector();

    void print();
  private:
    VectorType *Vector_;
    std::string filename_;
    /*std::vector<Vec3> C_;
    std::vector<Vec3> V_;
    int frame_;
    AtomMask mask_;
    AtomMask mask2_;*/

    int init();
    int setup();
    int action();
};
#endif
