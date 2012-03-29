#ifndef INC_ACTION_DISTRMSD_H
#define INC_ACTION_DISTRMSD_H
#include "Action.h"
#include "ActionReference.h"
// Class: DistRmsd
/// Action to calculate the distance RMSD between frame and a reference frame.
class DistRmsd: public Action, ActionReference {
  public:
    DistRmsd();

    int init();
    int setup();
    int action();
  private:
    DataSet *drmsd;    ///< DRMSD DataSet
    AtomMask TgtMask;  ///< Target mask.
    Frame SelectedTgt; ///< Hold only target coords selected by TgtMask
};
#endif
