#ifndef INC_ACTION_RMS2D_H
#define INC_ACTION_RMS2D_H
#include "Action.h"
#include "TrajectoryFile.h"
#include "CoordList.h"
#include "TriangleMatrix.h"
// Class: Action_Rms2d
/// Action to calculate the RMSD between two sets of frames.
/** Perform RMS calculation between each input frame and each other input 
  * frame, or each frame read in from a separate reference traj and each 
  * input frame. 
  * The actual calcuation is performed in the print function.
  */
class Action_Rms2d: public Action {
  public:
    Action_Rms2d();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Rms2d(); }
    static void Help();

    ~Action_Rms2d();

    void Print();
  private:
    CoordList ReferenceCoords_; ///< Hold coords from input frames.
    bool nofit_;                ///< Do not perform rms fitting
    bool useMass_;
    AtomMask RefMask_;          ///< Reference atom mask
    AtomMask FrameMask_;        ///< Target atom mask
    TrajectoryFile* RefTraj_;   ///< Reference trajectory, each frame used in turn
    Topology* RefParm_;         ///< Reference trajectory Parm
    Topology* mass_ptr_;        ///< If useMass, hold mass info for parm.
    DataSet* rmsdataset_;
    DataSet* Ct_;

    int AutoCorrelate(TriangleMatrix&);
    int CalcRmsToTraj();
    int Calc2drms();

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
};
#endif
