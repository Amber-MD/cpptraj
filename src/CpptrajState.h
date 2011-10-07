#ifndef INC_CPPTRAJSTATE_H
#define INC_CPPTRAJSTATE_H 
/// Class: CpptrajState
/// This is the main class for cpptraj. It holds all data and controls the 
/// overall flow of the program. It exists in main.cpp.
#include "TrajinList.h"
#include "TrajoutList.h"
#include "ReferenceList.h"
#include "ParmFileList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ActionList.h"
#include "AnalysisList.h"
#include "ArgList.h"
class CpptrajState {
    TrajinList trajinList;       // List of input trajectory files 
    ReferenceList referenceList; // List of reference coordinate files
    TrajoutList trajoutList;     // List of output trajectory files 
    ActionList actionList;       // List of actions to be performed each frame
    AnalysisList analysisList;   // List of analyses to be performed on datasets
    DataSetList DSL;             // List of generated data sets
    DataFileList DFL;            // List of datafiles that data sets will be written to
    std::list<ArgList*> DF_Args; // List of commands pertaining to datafile creation etc
    int TotalErrors;             // Sum of all returned error statuses
    int debug;

    void ProcessDataFileCmd();   // Handle datafile commands in DF_Args
    bool showProgress;           // Output traj progress to screen?

  public:
    void SetGlobalDebug(int);    // Set debug level for all components
    ParmFileList parmFileList;   // List of parameter files 
    void Dispatch(char*);        // Function that decides where to send commands

    CpptrajState();
    ~CpptrajState();

    int Run(); 
};
#endif
