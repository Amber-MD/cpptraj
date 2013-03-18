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

    void Print();
  private:
    enum vectorMode {
      NO_OP=0,   PRINCIPAL_X, PRINCIPAL_Y, PRINCIPAL_Z,
      DIPOLE,    BOX,         MASK,        IRED,
      CORRPLANE, CENTER 
    };
    static const char* ModeString[];

    DataSet_Vector* Vec_;
    DataSet* Magnitude_;
    double* vcorr_;
    vectorMode mode_;
    bool ptrajoutput_;
    Topology* CurrentParm_;
    AtomMask mask_;
    AtomMask mask2_;
    std::string filename_;

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    static double solve_cubic_eq(double,double,double,double);
    static Vec3 leastSquaresPlane(int,const double*);

    void Mask(Frame const&);
    void Dipole(Frame const&);
    void Principal(Frame const&);
    void CorrPlane(Frame const&);
};
#endif
