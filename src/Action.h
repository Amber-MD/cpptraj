#ifndef INC_ACTION_H
#define INC_ACTION_H
/* Action is the base class that all other actions inherit.
 */
#include "FrameList.h"
#include "ArgList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "AmberParm.h"

//enum actionType {
//  UNKNOWN_ACTION, ATOMMAP, DIHEDRAL, RMSD, DISTANCE, ANGLE, STRIP
//};

class Action {
  protected:
    ArgList *A;             // The action arguments (setArg)
    AmberParm *P;           // The current parmtop (setup)
    Frame *F;               // The current frame (action)
    DataSetList *DSL;       // Pointer to the data set list in PtrajState (init)
    DataFileList *DFL;      // Pointer to the data file list in PtrajState (init)
    FrameList *FL;          // Pointer to the reference frame list in PtrajState (init)
    int debug;
    int currentFrame;       // current frame being processed, set by PtrajActionList
    // --== Inherited by child classes ==--
    virtual int setup()  { return 0; }
    virtual int action() { return 0; }
    virtual int init() { return 0; }

  public:
    //actionType currentType; // The action type
    int noInit;             // Set to 1 if action could not be initialized
    int noSetup;            // Set to 1 if action could not be set up
  
    Action();               // Constructor
    virtual ~Action();      // Destructor - virtual since this class is inherited

    void setArg(ArgList *inA); // Set the argument list
    void ResetArg();           // Reset arguments in the argument list
    char *Name();              // Print the command that calls the action
    char *CmdLine();           // Print the entire argument line

    int Init(DataSetList*, FrameList*, DataFileList*, int); 
    int Setup(AmberParm **);
    int DoAction(Frame **,int);    

    // --== Inherited by child classes ==--
    virtual void info() {return;}
    virtual void print() { }   // Print anything besides datasets, called at end of execution 
};
#endif  
