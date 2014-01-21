#ifndef INC_ACTION_VECTOR_H
#define INC_ACTION_VECTOR_H
#include "Action.h"
#include "DataSet_Vector.h"
#include "TrajectoryFile.h"
class Action_Vector : public Action {
  public:
    Action_Vector();
    ~Action_Vector();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Vector(); }
    static void Help();
  private:
    enum vectorMode {
      NO_OP=0,   PRINCIPAL_X, PRINCIPAL_Y, PRINCIPAL_Z,
      DIPOLE,    BOX,         MASK,        IRED,
      CORRPLANE, CENTER 
    };
    static const char* ModeString[];

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    static double solve_cubic_eq(double,double,double,double);
    static Vec3 leastSquaresPlane(int,const double*);
    void Mask(Frame const&);
    void Dipole(Frame const&);
    void Principal(Frame const&);
    void CorrPlane(Frame const&);

    int ensembleNum_;
    DataSet_Vector* Vec_;   ///< Hold vector values
    DataSet* Magnitude_;    ///< Hold vector magnitudes if requested
    double* vcorr_;         ///< Temp. space for calculating CorrPlane
    vectorMode mode_;       ///< Vector calculation mode
    bool ptrajoutput_;      ///< If true output in ptraj format 
    Topology* CurrentParm_; ///< Current topology (for dipole)
    AtomMask mask_;
    AtomMask mask2_;
    std::string filename_;
    std::string trajoutName_;
    std::string parmoutName_;
    TrajectoryFile::TrajFormatType trajoutFmt_;
};
#endif
