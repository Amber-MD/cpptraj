#ifndef INC_CPPTRAJSTATE_H
#define INC_CPPTRAJSTATE_H
#include "TrajinList.h"
#include "TrajoutList.h"
#include "FrameList.h"
#include "TopologyList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ActionList.h"
#include "AnalysisList.h"
/// Hold all cpptraj state data
class CpptrajState {
  public:
    CpptrajState() : debug_(0), showProgress_(true), nrun_(0) {}
    // TODO: Change to &
    TopologyList* PFL()     { return &parmFileList_; }
    FrameList* FL()         { return &refFrames_;    }
    DataSetList* DSL()      { return &DSL_;          }
    DataFileList* DFL()     { return &DFL_;          }
    ActionList& ActList()   { return actionList_;    }
    AnalysisList& AnaList() { return analysisList_;  }
    int Debug()       const { return debug_;         }
    //ArgsArray& TrajoutArgs() { return trajoutArgs_; }
    //ArgsArray& ActionArgs()  { return actionArgs_; }
    int ListAll(ArgList&);
    int SetListDebug(ArgList&);
    int ClearList(ArgList&);
    int Run();
  private:
    /// Types of lists
    enum ListType {
      L_ACTION = 0, L_TRAJIN, L_REF, L_TRAJOUT, L_PARM, L_ANALYSIS,
      L_DATAFILE, L_DATASET, N_LISTS
    };
    std::vector<bool> ListsFromArg(ArgList&, bool);

    int RunNormal();
    int RunEnsemble();

    /// List of parameter files 
    TopologyList parmFileList_;
    /// List of input trajectory files
    TrajinList trajinList_;
    /// List of reference coordinate files
    FrameList refFrames_;
    
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
    
    typedef std::vector<ArgList> ArgsArray;
    /// Array of trajout args for setting up ensemble trajout.
    ArgsArray trajoutArgs_;
    /// Array of action args for setting up ensemble actions.
    //ArgsArray actionArgs_;
    /// State debug level
    int debug_;
    /// Display Progress bar during run
    bool showProgress_;
    /// Number of times Run() has been called.
    int nrun_;
};
#endif
