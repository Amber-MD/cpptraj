#ifndef INC_CPPTRAJSTATE_H
#define INC_CPPTRAJSTATE_H 
#include "TrajinList.h"
#include "TrajoutList.h"
#include "ReferenceList.h"
#include "ParmFileList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ActionList.h"
#include "AnalysisList.h"
#include "ArgList.h"
// Class: CpptrajState:
/// Hold state information.
/** This is the main class for cpptraj. It holds all data and controls the 
 *  overall flow of the program. It exists in main.cpp.
 */
class CpptrajState {
    /// List of input trajectory files
    TrajinList trajinList;
    /// List of reference coordinate files
    ReferenceList referenceList; 
    /// List of output trajectory files 
    TrajoutList trajoutList;
    /// List of actions to be performed each frame
    ActionList actionList;    
    /// List of analyses to be performed on datasets
    AnalysisList analysisList;
    /// List of generated data sets
    DataSetList DSL;
    /// List of datafiles that data sets will be written to
    DataFileList DFL;
    /// The debug level
    int debug;
    /// The number of the active reference structure in ReferenceList
    int activeRef;
    /// If true the progress of reading input trajectories will be shown
    bool showProgress;
    /// If true cpptraj will exit if errors are encountered instead of trying to continue
    bool exitOnError;
  public:
    /// Set debug level for all components
    void SetGlobalDebug(int);
    /// List of parameter files 
    ParmFileList parmFileList;
    /// Function that decides where to send commands
    void Dispatch(char*);        

    CpptrajState();
    ~CpptrajState();
    /// Controls main flow of the program.
    int Run(); 
};
#endif
