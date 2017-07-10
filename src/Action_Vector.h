#ifndef INC_ACTION_VECTOR_H
#define INC_ACTION_VECTOR_H
#include "Action.h"
#include "DataSet_Vector.h"
class Action_Vector : public Action {
  public:
    Action_Vector();
    ~Action_Vector();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Vector(); }
    void Help() const;
  private:
    enum vectorMode {
      NO_OP=0,   PRINCIPAL_X, PRINCIPAL_Y, PRINCIPAL_Z,
      DIPOLE,    BOX,         MASK,        IRED,
      CORRPLANE, CENTER,      BOX_X,       BOX_Y,       BOX_Z,
      BOX_CTR,   MINIMAGE,    MOMENTUM
    };
    static const char* ModeString[];

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    static double solve_cubic_eq(double,double,double,double);
    static Vec3 leastSquaresPlane(int,const double*);
    void Mask(Frame const&);
    void Dipole(Frame const&);
    void Principal(Frame const&);
    void CorrPlane(Frame const&);
    void UnitCell(Box const&);
    void MinImage(Frame const&);

    DataSet_Vector* Vec_;   ///< Hold vector values
    DataSet* Magnitude_;    ///< Hold vector magnitudes if requested
    double* vcorr_;         ///< Temp. space for calculating CorrPlane
    vectorMode mode_;       ///< Vector calculation mode
    bool ptrajoutput_;      ///< If true output in ptraj format
    bool needBoxInfo_;      ///< If true box info required. 
    Topology* CurrentParm_; ///< Current topology (for dipole)
    AtomMask mask_;
    AtomMask mask2_;
    CpptrajFile* outfile_;
};
#endif
