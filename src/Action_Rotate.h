#ifndef INC_ACTION_ROTATE_H
#define INC_ACTION_ROTATE_H
#include "Action.h"
#include "DataSet_Mat3x3.h"
class Action_Rotate : public Action {
  public:
    Action_Rotate();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Rotate(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    enum ModeType { ROTATE = 0, DATASET, AXIS };
    Matrix_3x3 RotMatrix_;      ///< Rotation matrix.
    AtomMask mask_;             ///< Mask of atoms to rotate.
    AtomMask axis0_;            ///< Mask of atoms defining 1 end of rotation axis.
    AtomMask axis1_;            ///< Mask of atoms defining other end of rotation axis.
    DataSet_Mat3x3* rmatrices_; ///< DataSet containing rotation matrices to use.
    double delta_;              ///< Degrees to rotate if mode is AXIS
    ModeType mode_;             ///< Mode to use.
    bool inverse_;              ///< If true perform an inverse rotation.
};
#endif
