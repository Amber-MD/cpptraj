#ifndef INC_ACTION_VECTOR_H
#define INC_ACTION_VECTOR_H
#include "Action.h"
#include "DataSet_Vector.h"
class Action_Vector : public Action {
  public:
    Action_Vector();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Vector(); }
    static void Help();

    ~Action_Vector();

    void print();
  private:
    enum vectorMode {
      NO_OP=0,   PRINCIPAL_X, PRINCIPAL_Y, PRINCIPAL_Z,
      DIPOLE,    BOX,         MASK,        IRED,
      CORRPLANE, CORR,        CORRIRED
    };
    static const char ModeString[][12];

    DataSet_Vector* Vec_;
    double* vcorr_;
    vectorMode mode_;
    bool ptrajoutput_;
    AtomMask mask_;
    AtomMask mask2_;
    std::string filename_;

    int init();
    int setup();
    int action();

    static double solve_cubic_eq(double,double,double,double);
    static Vec3 leastSquaresPlane(int,const double*);

    void Mask();
    void Dipole();
    void Principal();
    void CorrPlane();
    void Box();
};
#endif
