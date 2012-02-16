#ifndef INC_ACTION_H
#define INC_ACTION_H
// Class: Action 
/// The base class that all other actions inherit. 
/** Each action currently
  * holds the memory address of all important state information: the
  * master dataset list, datafile list, parm file list, and reference
  * frame list. By convention actions have 4 main phases: Init, Setup,
  * Action, and print (optional).
  * Mass Common Functionality:
  * Since several actions have the option to include mass information as
  * part of the calculation (mostly for calculating center of mass as opposed 
  * to geometric center) the variable useMass is common to all actions and
  * mass information is checked in Setup. useMass is set to false if the parm
  * does not contain mass information.
  */
// NOTE: This means that mass disabled for all future invocations of setup,
//       should probably be changed.
#include "FrameList.h"
#include "ArgList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ParmFileList.h"
#include "AmberParm.h"
class Action {
  protected:
    bool isSeparate;        ///< If true action was initialized outside main action list.
    ArgList actionArgs;     ///< The action arguments (setArg)
    AmberParm *currentParm; ///< The current parmtop (setup)
    Frame *currentFrame;    ///< The current frame (action)
    DataSetList *DSL;       ///< Pointer to the data set list in CpptrajState (init)
    DataFileList *DFL;      ///< Pointer to the data file list in CpptrajState (init)
    ParmFileList *PFL;      ///< Pointer to the parm file list in CpptrajState (init)
    FrameList *FL;          ///< Pointer to the reference frame list in CpptrajState (init)
    Frame *activeRefFrame;  ///< Pointer to active ref frame in frame list
    double *activeReference; ///< Temp pointer to coords of active ref frame in frame list
    /// Pointer to coords of active ref frame in frame list
    double *activeReferenceOriginalValue;
    
    bool useMass;           ///< If set to true, calculations will use mass info
    bool useMassOriginalValue; ///< Value of useMass set by init

    int debug;              ///< action debug level
    int frameNum;           ///< # of current frame being processed, set by ActionList
    // --== Inherited by child classes ==--
    /// actions internal setup routine, called by Setup
    /** Setup action. Process any parm-dependent things like masks.
      */
    virtual int setup()  { return 0; }
    /// actions internal action routine, called by DoAction
    /** Perform action. Only parts of the action which depend on input
      * coordinates should be implemented here.
      */
    virtual int action() { return 0; }
    /// actions internal init routine, called by Init
    /** Initialize action. Parse arguments from actionArgs. Called prior to reading
      * input trajectories. Expected order of checking args is: 1) Keywords, 2) 
      * Masks, 3) dataset names.
      */
    virtual int init()   { return 0; }

  public:
    ///< Enumerate potential return states from DoAction.
    enum ActionReturnType { ACTION_OK=0, ACTION_ERR, ACTION_USEORIGINALFRAME,
                            ACTION_SUPPRESSCOORDOUTPUT 
                          };
    bool noInit;             ///< True if action could not be initialized
    bool noSetup;            ///< True if action could not be set up
  
    Action();               // Constructor
    virtual ~Action();      // Destructor - virtual since this class is inherited

    void SetArg(const ArgList&); ///< Set the argument list
    const char *ActionCommand(); ///< Print the command that calls the action
    const char *CmdLine();       ///< Print the entire argument line

    /// Initialize action, call init()
    int Init(DataSetList*, FrameList*, DataFileList*, ParmFileList *,int); 
    /// Set up action for parm, call setup()
    int Setup(AmberParm **);
    /// Perform action on frame, call action()
    ActionReturnType DoAction(Frame **,int);    

    // --== Inherited by child classes ==--
    /// Print anything besides datasets, called at end of execution
    /** Perform any output not related to master dataset output, or any 
      * necessary post-trajectory calculations, e.g. in the 2drms command
      * several passes must be made over the input frames, which are stored
      * by that action.
      */
    virtual void print() { } 
};
#endif  
