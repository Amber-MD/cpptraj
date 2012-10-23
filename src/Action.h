#ifndef INC_ACTION_H
#define INC_ACTION_H
#include "DispatchObject.h"
#include "ArgList.h"
#include "DataFileList.h"
#include "DataSetList.h"
#include "FrameList.h"
#include "TopologyList.h"
// Class: Action 
/// The abstract base class that all other actions inherit. 
/** By convention actions have 3 main phases: Init, Setup, and DoAction.
  * Init is used to initialize the action, make sure that all arguments
  * for the action are correct, and add any DataSets/DataFiles which will
  * be used by the action. Setup will set up the action for a specific
  * Topology file. DoAction will perform the action on a given frame.
  * A fourth function, Print, is for any additional calculations or output 
  * the action may require once all frames are processed.
  */
class Action : public DispatchObject {
  public:
    /// Enumerate potential return states from Init, Setup, and DoAction.
    enum RetType { OK=0, ERR, USEORIGINALFRAME, SUPPRESSCOORDOUTPUT };
    /// Destructor - virtual since this class is inherited
    virtual ~Action() {}
    /// Initialize action
    /** Process input args, set up any DataSets or DataFiles, set debug level */
    virtual RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*, 
                         DataFileList*, int) = 0;
    /// Set up action for given Topology
    virtual RetType Setup(Topology*,Topology**) = 0;
    /// Perform action for given frame number and Frame.
    virtual RetType DoAction(int,Frame*,Frame**) = 0;
    /// Print anything besides datasets, called at end of execution
    /** Perform any output not related to master dataset output, or any 
      * necessary post-trajectory calculations, e.g. in the 2drms command
      * several passes must be made over the input frames, which are stored
      * by that action.
      */
    virtual void Print() = 0;
};
#endif  
