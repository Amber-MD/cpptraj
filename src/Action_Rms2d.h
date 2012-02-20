#ifndef INC_ACTION_RMS2D_H
#define INC_ACTION_RMS2D_H
#include "Action.h"
#include "TrajectoryFile.h"
#include "CoordList.h"
#include "TriangleMatrix.h"
#include "DataSet_double.h"
// Class: Rms2d
/// Action to calculate the RMSD between two sets of frames.
/** Perform RMS calculation between each input frame and each other input 
  * frame, or each frame read in from a separate reference traj and each 
  * input frame. 
  * The actual calcuation is performed in the print function.
  */
class Rms2d: public Action {
    CoordList ReferenceCoords; ///< Hold coords from input frames.
    bool nofit;                ///< Do not perform rms fitting
    AtomMask RefMask;          ///< Reference atom mask
    AtomMask FrameMask;        ///< Target atom mask
    char *rmsdFile;            ///< Output filename
    DataSetList RmsData;       ///< 1 data set for each ref frame to each tgt frame
    TrajectoryFile *RefTraj;   ///< Reference trajectory, each frame used in turn
    AmberParm *RefParm;        ///< Reference trajectory Parm
    DataSet_double Ct;         ///< Hold auto-correlation
    char *corrfilename;        ///< Auto-correlation output filename
    double *mass_ptr;          ///< If useMass, hold mass info for parm.
    bool mass_setup;           ///< Used to check if mass already set up.

    void CalcRmsToTraj();
    int AutoCorrelate(TriangleMatrix &);
  public:
    Rms2d();
    ~Rms2d();

    int SeparateInit(bool, char *);
    void Calc2drms(TriangleMatrix*);

    int init();
    int setup();
    int action();
    void print();
};
#endif
