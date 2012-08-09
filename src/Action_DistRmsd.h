#ifndef INC_ACTION_DISTRMSD_H
#define INC_ACTION_DISTRMSD_H
#include "Action.h"
#include "TrajectoryFile.h"
// Class: Action_DistRmsd
/// Action to calculate the distance RMSD between frame and a reference frame.
class Action_DistRmsd: public Action {
  public:
    Action_DistRmsd();
  private:
    int init();
    int setup();
    int action();

    DataSet *drmsd_;    ///< DRMSD DataSet
    AtomMask TgtMask_;  ///< Target mask.
    Frame SelectedTgt_; ///< Hold only target coords selected by TgtMask
    // Reference variables and functions
    enum RefModeType { UNKNOWN_REF=0, FIRST, REF, REFTRAJ };
    RefModeType refmode_;
    Frame RefFrame_;
    Frame SelectedRef_;
    AtomMask RefMask_;
    TrajectoryFile RefTraj_;

    int SetRefMask( Topology* );
    void SetRefStructure( Frame& );
};
#endif
