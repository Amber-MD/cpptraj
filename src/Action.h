#ifndef INC_ACTION_H
#define INC_ACTION_H
/// Class: Action 
/// The base class that all other actions inherit. Each action currently
/// holds the memory address of all important state information: the
/// master dataset list, datafile list, parm file list, and reference
/// frame list. By convention actions have 4 main phases:
///   Init: 
///     Initialize action, get all relevant arguments from the
///     argument list, set up any datasets which will be available for 
///     analysis after trajectory processing, and set memory references 
///     to the master dataset, datafile, reference frame, and parm file lists.
///     Actions can have datasets not related to master datasets, e.g.
///     in the NAstruct action, since data is stored for each nucleobase at
///     each frame its output would not match up with other actions in the 
///     master dataset list, so it has its own dataset list.
///   Setup:
///     Set up action for the current parm file. This is where any 
///     parm-dependent variables should be set such as atom masks etc. The
///     current parm memory address is set here but can also be modified
///     by the action, this allows e.g. stripping of the parm. Only copies
///     of the parm should be modified; a reference to the original parm is
///     always stored in CpptrajState and can be reset there with the 'unstrip'
///     command.
///   Action:
///     Perform action on the current frame. The current frame memory address
///     is passed in and can be modified by the action, again for things like
///     stripping etc. The current frame number is also passed in.
///   Print (optional):
///     Perform any output not related to master dataset output, or any 
///     necessary post-trajectory calculations, e.g. in the 2drms command
///     several passes must be made over the input frames, which are stored
///     by that action.
/// Mass Common Functionality:
/// Since several actions have the option to include mass information as
/// part of the calculation (mostly for calculating center of mass as opposed 
/// to geometric center) the variable useMass is common to all actions and
/// mass information is checked in Setup. useMass is set to false if the parm
/// does not contain mass information.
/// NOTE: This means that mass disabled for all future invocations of setup,
///       should probably be changed.
#include "FrameList.h"
#include "ArgList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ParmFileList.h"
#include "AmberParm.h"
class Action {
  protected:
    ArgList actionArgs;     // The action arguments (setArg)
    AmberParm *currentParm; // The current parmtop (setup)
    Frame *currentFrame;    // The current frame (action)
    DataSetList *DSL;       // Pointer to the data set list in CpptrajState (init)
    DataFileList *DFL;      // Pointer to the data file list in CpptrajState (init)
    ParmFileList *PFL;      // Pointer to the parm file list in CpptrajState (init)
    FrameList *FL;          // Pointer to the reference frame list in CpptrajState (init)
    double *activeReference; // Pointer to coords of active ref frame in frame list
    
    bool useMass;           // If set to true, calculations will use mass info

    int debug;
    int frameNum;           // # of current frame being processed, set by ActionList
    // --== Inherited by child classes ==--
    virtual int setup()  { return 0; }
    virtual int action() { return 0; }
    virtual int init()   { return 0; }

  public:
    int noInit;             // Set to 1 if action could not be initialized
    int noSetup;            // Set to 1 if action could not be set up
  
    Action();               // Constructor
    virtual ~Action();      // Destructor - virtual since this class is inherited

    void SetArg(const ArgList&); // Set the argument list
    const char *ActionCommand();              // Print the command that calls the action
    const char *CmdLine();           // Print the entire argument line

    int Init(DataSetList*, FrameList*, DataFileList*, ParmFileList *,int); 
    int Setup(AmberParm **);
    int DoAction(Frame **,int);    

    // --== Inherited by child classes ==--
    virtual void print() { }   // Print anything besides datasets, called at end of execution 
};
#endif  
