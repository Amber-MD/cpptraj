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
// Class: Cpptraj:
/// Hold state information.
/** This is the main class for cpptraj. It holds all data and controls the 
 *  overall flow of the program. It exists in main.cpp.
 */
class Cpptraj {
  public:
    enum Mode { C_OK = 0, C_ERR, C_QUIT, C_INTERACTIVE };
    Cpptraj();
    Mode Interactive();
    Mode ProcessInput(std::string const&);
    Mode ProcessCmdLineArgs(int,char**);
    int Run();
  private:
    static void Usage();

    void Help(ArgList&);
    void ListAction(ArgList&,int);
    int Create_DataFile(ArgList&,int);
    int Precision(ArgList&);
    int ReadData(ArgList&);
    void SelectDS(ArgList&);
    int LoadParm(ArgList&);
    int ParmInfo(ArgList&,int);
    int ParmWrite(ArgList&);
    int ParmStrip(ArgList&);
    int ParmBox(ArgList&);
    int ParmSolvent(ArgList&);
    int Select(ArgList&);
    int LoadCrd(ArgList&);
    int CrdAction(ArgList&);
    int CrdOut(ArgList&);
    int CrdAnalyze(ArgList&);
    /// Function that decides where to send commands
    Mode Dispatch(std::string const&);
    /// List of parameter files 
    TopologyList parmFileList_;
    /// List of input trajectory files
    TrajinList trajinList_;
    /// List of reference coordinate files
    FrameList refFrames_; 

    typedef std::vector<ArgList> ArgsArray;
    /// Array of trajout args for setting up ensemble trajout.
    ArgsArray trajoutArgs_;
    /// Array of action args for setting up ensemble actions.
    ArgsArray actionArgs_;

    /// List of output trajectory files 
    TrajoutList trajoutList_;
    /// List of analyses to be performed on datasets
    AnalysisList analysisList_;
    /// List of generated data sets
    DataSetList DSL_;
    /// List of actions to be performed each frame
    ActionList actionList_;    
    /// List of datafiles that data sets will be written to
    DataFileList DFL_;
    /// The debug level
    int debug_;
    /// If true the progress of reading input trajectories will be shown
    bool showProgress_;
    /// If true cpptraj will exit if errors are encountered instead of trying to continue
    bool exitOnError_;
    /// Number of times the Run routine has been called.
    int nrun_;
    /// Log file for interactive mode
    CpptrajFile logfile_;

    int RunEnsemble();
    int RunNormal();
};
#endif
