#ifndef INC_CPPTRAJSTATE_H
#define INC_CPPTRAJSTATE_H
#include "TrajinList.h"
#include "TrajoutList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ActionList.h"
#include "AnalysisList.h"
/// Hold all cpptraj state data
class CpptrajState {
  public:
    /// Possible command return types. Put here because both Command and Exec need this.
    enum RetType { OK = 0, ERR, QUIT };
    /// CONSTRUCTOR
    CpptrajState() : debug_(0), showProgress_(true), exitOnError_(true), noEmptyRun_(false) {}
    void SetNoExitOnError()  { exitOnError_ = false;  }
    void SetNoProgress()     { showProgress_ = false; }
    void SetActionSilence(bool b)  { actionList_.SetSilent(b); }

    DataSetList const& DSL()  const { return DSL_;         }
    DataSetList&       DSL()        { return DSL_;         }
    DataFileList const& DFL() const { return DFL_;         }
    DataFileList&       DFL()       { return DFL_;         }
    int Debug()               const { return debug_;       }
    bool ExitOnError()        const { return exitOnError_; }
    bool EmptyState()         const { return (actionList_.Empty() && 
                                              analysisList_.Empty() &&
                                              trajoutList_.Empty()); }
    TrajinList const& InputTrajList() const { return trajinList_; }

    int AddInputTrajectory( ArgList&, bool );
    int AddInputTrajectory( std::string const& );
    int AddOutputTrajectory( ArgList& );
    int AddOutputTrajectory( std::string const& );
    int RunAnalyses();
    // TODO: Move AddReference() to DataSetList?
    int AddReference( std::string const&, ArgList const& );
    inline int AddReference( std::string const& );
    int AddTopology( std::string const&, ArgList const& );
    int AddTopology( Topology const&, std::string const& );
    inline int AddToActionQueue( Action*, ArgList& );
    inline int AddToAnalysisQueue( Analysis*, ArgList& );
    static int WorldSize();
    static std::string PrintListKeys();
    int ListAll(ArgList&) const;
    int SetListDebug(ArgList&);
    int ClearList(ArgList&);
    int RemoveDataSet(ArgList&);
    int TrajLength( std::string const&, std::vector<std::string> const&);
    int Run();
    /// Write all DataFiles
    void MasterDataFileWrite();
  private:
    /// Types of lists
    enum ListType {
      L_ACTION = 0, L_TRAJIN, L_REF, L_TRAJOUT, L_PARM, L_ANALYSIS,
      L_DATAFILE, L_DATASET, N_LISTS
    };
    /// Hold list keyword.
    struct ListKeyType {
      ListType Type_;
      const char* Key_;
    };
    static ListKeyType ListKeys[];
    std::vector<bool> ListsFromArg(ArgList&, bool) const;

    int RunNormal();
    int RunEnsemble();
#   ifdef MPI
    std::vector<int> DivideFramesAmongThreads(int&, int&, int&, int, int, int, bool);
    int RunParallel();
    //int RunSingleTrajParallel();
#   endif
    // -------------------------------------------
     /// List of generated data sets
    DataSetList DSL_;
    /// List of datafiles that data sets will be written to
    DataFileList DFL_;
    /// List of input trajectory files
    TrajinList trajinList_;
    // -------------------------------------------
    /// List of actions to be performed each frame
    ActionList actionList_;
    /// List of output trajectory files 
    TrajoutList trajoutList_;
    // -------------------------------------------
    /// List of analyses to be performed on datasets
    AnalysisList analysisList_;
    
    /// State debug level
    int debug_;
    /// Display Progress bar during run
    bool showProgress_;
    /// If true cpptraj will exit if errors are encountered instead of trying to continue
    bool exitOnError_;
    /// If true do not process input trajectories when no actions/output trajectories.
    bool noEmptyRun_; // DEBUG: false is used for benchmarking trajectory read speed.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// CpptrajState::AddToActionQueue()
int CpptrajState::AddToActionQueue( Action* actIn, ArgList& argIn ) {
  argIn.MarkArg(0);
  ActionInit init(DSL_, DFL_);
  return actionList_.AddAction( actIn, argIn, init );
}
// CpptrajState::AddToAnalysisQueue()
int CpptrajState::AddToAnalysisQueue( Analysis* anaIn, ArgList& argIn ) {
  argIn.MarkArg(0);
  AnalysisSetup setup(DSL_, DFL_);
  return analysisList_.AddAnalysis( anaIn, argIn, setup );
}
// CpptrajState::AddReference()
int CpptrajState::AddReference( std::string const& fname ) {
  return AddReference( fname, ArgList() );
}
#endif
