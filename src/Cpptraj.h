#ifndef INC_CPPTRAJ_H
#define INC_CPPTRAJ_H 
#include "TrajinList.h"
#include "TrajoutList.h"
#include "FrameList.h"
#include "TopologyList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ActionList.h"
#include "AnalysisList.h"
#include "ArgList.h"
// Class: Cpptraj:
/// Hold state information.
/** This is the main class for cpptraj. It holds all data and controls the 
 *  overall flow of the program. It exists in main.cpp.
 */
class Cpptraj {
  public:
    /// Set debug level for all components
    void SetGlobalDebug(int);
    void AddParm(const char *);
    void Interactive();
    /// Function that decides where to send commands
    bool Dispatch(const char*);        

    Cpptraj();
    /// Controls main flow of the program.
    int Run();
  private:
    static void Help_List();
    static void Help_Help();
    static void Help_Debug();
    static void Help_ActiveRef();
    void List(ArgList&);
    void Help(ArgList&);
    void Debug(ArgList&);

    static const DispatchObject::Token GeneralCmds[];
    static const DispatchObject::Token CoordCmds[];
    int SearchTokenArray(const DispatchObject::Token* DispatchArray,
                         bool, const ArgList&);
    int SearchToken(const ArgList&);
    DispatchObject::Token const* dispatchToken_;
    /// List of parameter files 
    TopologyList parmFileList;
    /// List of input trajectory files
    TrajinList trajinList;
    /// List of reference coordinate files
    FrameList refFrames; 
    /// List of output trajectory files 
    TrajoutList trajoutList;
    /// List of analyses to be performed on datasets
    AnalysisList analysisList;
    /// List of generated data sets
    DataSetList DSL;
    /// List of actions to be performed each frame
    ActionList actionList;    
    /// List of datafiles that data sets will be written to
    DataFileList DFL;
    /// The debug level
    int debug_;
    /// If true the progress of reading input trajectories will be shown
    bool showProgress_;
    /// If true cpptraj will exit if errors are encountered instead of trying to continue
    bool exitOnError_;
};
#endif
