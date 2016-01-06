#ifndef INC_ACTION_H
#define INC_ACTION_H
#include "DispatchObject.h"
#include "ArgList.h"
#include "ActionState.h"
#ifdef MPI
# include "Parallel.h"
#endif
/// The abstract base class that all other Actions inherit.
/** By convention Actions have 3 main phases: Init, Setup, and DoAction.
  * Init is used to initialize the Action, make sure that all arguments
  * for the Action are correct, and add any DataSets/DataFiles which will
  * be used by the Action. Setup will set up the Action for a specific
  * Topology file. DoAction will perform the Action on a given frame.
  * A fourth function, Print, is for any additional calculations or output 
  * the Action may require once all frames are processed.
  * Note that Setup and DoAction are allowed to modify the passed in Topology
  * and Frame - they do so by creating their own copies (and in such cases they
  * must return MODIFY_TOPOLOGY/MODIFY_COORDS respectively). This allows
  * Actions to modify the current Topology/Frame by replacing them with ones
  * inside the Action itself. This is more memory-hungry but allows the
  * Topology/Frame to be manipulated quickly, and is also easier to return to
  * the original Topology/Frame via USE_ORIGINAL_FRAME (see e.g. Action_Unstrip
  * in Action_Unstrip.cpp).
  */
class Action : public DispatchObject {
  public:
    /// Constructor
    Action() : DispatchObject(ACTION) {}
    /// Constructor - override ACTION (e.g. HIDDEN)
    Action(DispatchObject::Otype o) : DispatchObject(o) {}
    /// Enumerate potential return states from Init, Setup, and DoAction.
    enum RetType { OK=0, ///< Everything OK, normal return.
                   ERR,  ///< Problem occurred.
                   USE_ORIGINAL_FRAME, ///< Return to unmodified frame/topology.
                   SUPPRESS_COORD_OUTPUT, ///< Skip remaining Actions and traj output. TODO just SKIP?
                   SKIP, ///< Non-fatal problem occurred, skip Action until re-setup.
                   MODIFY_TOPOLOGY, ///< Action has modified the topology.
                   MODIFY_COORDS ///< Action has modified the frame.
    };
    /// Destructor - virtual since this class is inherited
    virtual ~Action() {}
    /// Initialize Action
    /** Process input args, set up any DataSets or DataFiles, set debug level */
    virtual RetType Init(ArgList&, ActionInit&, int) = 0;
    /// Set up Action for given Topology
    virtual RetType Setup(ActionSetup&) = 0;
    /// Perform Action for given frame number and Frame.
    virtual RetType DoAction(int, ActionFrame&) = 0;
    /// Print anything besides datasets, called at end of execution
    /** Perform any output not related to master dataset output, or any 
      * necessary post-trajectory processing calculations.
      */
    virtual void Print() = 0;
#   ifdef MPI
    /// Sync Action data to master when running in parallel across trajectories.
    virtual int SyncAction(Parallel::Comm const&) { return 0; } // TODO: pure virtual
#   endif
};
#endif
