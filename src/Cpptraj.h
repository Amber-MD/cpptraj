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
    Mode ProcessCmdLineArgs(int,char**);
    int Run();
  private:
    static void Usage(const char*);
    static void Help_Help();
    static void Help_List();
    static void Help_Debug();
    static void Help_Clear();
    static void Help_ActiveRef();
    static void Help_Create_DataFile();
    static void Help_Precision();
    static void Help_SelectDS();
    static void Help_Trajin();
    static void Help_Ensemble();
    static void Help_Trajout();
    static void Help_Reference();
    static void Help_Parm();
    static void Help_ParmInfo();
    static void Help_ParmWrite();
    static void Help_ParmStrip();
    static void Help_ParmBox();
    static void Help_Solvent();
    static void Help_BondInfo();
    static void Help_ResInfo();
    static void Help_MolInfo();
    static void Help_CrdAction();
    static void Help_CrdOut();
    void Help(ArgList&);
    void List(ArgList&);
    void Debug(ArgList&);
    void Clear(ArgList&);
    int Create_DataFile(ArgList&);
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
    int CrdAction(ArgList&);
    int CrdOut(ArgList&);
    int CrdAnalyze(ArgList&);

    static const DispatchObject::Token GeneralCmds[];
    static const DispatchObject::Token CoordCmds[];
    static const DispatchObject::Token ParmCmds[];
    static const DispatchObject::Token Deprecated[];
    void ListAllCommands(const DispatchObject::Token*);
    DispatchObject::TokenPtr SearchTokenArray(const DispatchObject::Token*, ArgList const&);
    DispatchObject::TokenPtr SearchToken(ArgList&);

    Mode ProcessInput(std::string const&);
    Mode Dispatch(const char*);        ///< Function that decides where to send commands

    /// List of parameter files 
    TopologyList parmFileList;
    /// List of input trajectory files
    TrajinList trajinList;
    /// List of reference coordinate files
    FrameList refFrames; 

    typedef std::vector<ArgList> ArgsArray;
    /// Array of trajout args for setting up ensemble trajout.
    ArgsArray trajoutArgs_;
    /// Array of action args for setting up ensemble actions.
    ArgsArray actionArgs_;

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
    /// Number of times the Run routine has been called.
    int nrun_;
    /// Log file for interactive mode
    CpptrajFile logfile_;

    int RunEnsemble();
    int RunNormal();
};
#endif
