#ifndef INC_PTRAJSTATE_H
#define INC_PTRAJSTATE_H 

#include "TrajinList.h"
#include "TrajoutList.h"
#include "ReferenceList.h"
#include "ParmFileList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "PtrajActionList.h"
#include "ArgList.h"

class PtrajState {
    TrajinList trajFileList;         // List of input trajectory files 
    ReferenceList refFileList;       // List of reference coordinate files
    TrajoutList outFileList;         // List of output trajectory files 
    ParmFileList parmFileList;       // List of parameter files 
    PtrajActionList ptrajActionList; // List of actions to be performed each frame
    DataSetList DSL;                 // List of generated data sets
    DataFileList DFL;                // List of datafiles that data sets will be written to
    ArgList *A;                      // Current argument list
    std::list<ArgList*> DF_Args;     // List of commands pertaining to datafile creation etc
    int TotalErrors;                 // Sum of all returned error statuses
    int debug;

    void SetGlobalDebug(int);        // Set debug level for all components
    void Dispatch();                 // Function that decides where to send commands
    void ProcessDataFileCmd();       // Handle datafile commands in DF_Args
    int ProcessInputStream(char *);  // Process lines of input from a file
    int showProgress;                // Output traj progress to screen?

  public:

    PtrajState();
    ~PtrajState();

    int ProcessCmdLineArgs(int, char **); // Process arguments from command line
    int Run(); 
};
#endif
