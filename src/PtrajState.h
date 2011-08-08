#ifndef INC_PTRAJSTATE_H
#define INC_PTRAJSTATE_H 
/// Class: PtrajState
/// This is the main class for cpptraj. It holds all data and controls the 
/// overall flow of the program. It exists in main.cpp.
#include "TrajinList.h"
#include "TrajoutList.h"
#include "ReferenceList.h"
#include "ParmFileList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "PtrajActionList.h"
#include "AnalysisList.h"
#include "ArgList.h"
class PtrajState {
    TrajinList trajinList;           // List of input trajectory files 
    ReferenceList referenceList;     // List of reference coordinate files
    TrajoutList trajoutList;         // List of output trajectory files 
    ParmFileList parmFileList;       // List of parameter files 
    PtrajActionList ptrajActionList; // List of actions to be performed each frame
    AnalysisList analysisList;       // List of analyses to be performed on datasets
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
    bool showProgress;               // Output traj progress to screen?

  public:

    PtrajState();
    ~PtrajState();

    int ProcessCmdLineArgs(int, char **); // Process arguments from command line
    int Run(); 
};
#endif
