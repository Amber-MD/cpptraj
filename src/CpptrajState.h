#ifndef INC_CPPTRAJSTATE_H
#define INC_CPPTRAJSTATE_H
#include "TrajinList.h"
#include "TrajoutList.h"
#include "EnsembleOutList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ActionList.h"
#include "AnalysisList.h"
#include "Timer.h"
/// Hold all cpptraj state data
class CpptrajState {
  public:
    /// Possible command return types. Put here because both Command and Exec need this.
    enum RetType { OK = 0, ERR, QUIT };
    /// Trajectory mode
    enum TrajModeType { UNDEFINED = 0, NORMAL, ENSEMBLE };
    /// CONSTRUCTOR
    CpptrajState();
    void SetNoExitOnError()  { exitOnError_ = false;  }
    void SetNoProgress()     { showProgress_ = false; }
    void SetQuietBlocks(bool b)    { quietBlocks_ = b;         }
    void SetActionSilence(bool b)  { actionList_.SetSilent(b); }
#   ifdef MPI
    void SetForceParaEnsemble(bool b) { forceParallelEnsemble_ = b; }
#   endif
    DataSetList const& DSL()  const { return DSL_;         }
    DataSetList&       DSL()        { return DSL_;         }
    DataFileList const& DFL() const { return DFL_;         }
    DataFileList&       DFL()       { return DFL_;         }
    AnalysisList const& Analyses() const { return analysisList_; }
    AnalysisList&       Analyses()       { return analysisList_; }
    TrajModeType Mode()       const { return mode_;        }
    int Debug()               const { return debug_;       }
    bool ShowProgress()       const { return showProgress_;}
    bool QuietBlocks()        const { return quietBlocks_; }
    bool ExitOnError()        const { return exitOnError_; }
    bool RecordAllInput()     const { return recordAllInput_; }
    bool EmptyState()         const { return (actionList_.Empty() && 
                                              analysisList_.Empty() &&
                                              trajoutList_.Empty() &&
                                              ensembleOut_.Empty()); }
    TrajinList const& InputTrajList() const { return trajinList_; }

    int AddInputTrajectory( std::string const& );
    int AddInputTrajectory( ArgList& );
    int AddInputEnsemble( ArgList& );
    int AddOutputTrajectory( ArgList& );
    int AddOutputTrajectory( std::string const& );
    int RunAnalyses();
    // TODO: Move AddReference() to DataSetList?
    int AddReference( std::string const&, ArgList const& );
    int AddReference( std::string const& );
    int AddTopology( std::string const&, ArgList const& );
    int AddTopology( Topology const&, std::string const& );
    RetType AddToActionQueue( Action*, ArgList& );
    RetType AddToAnalysisQueue( Analysis*, ArgList& );
    static std::string PrintListKeys();
    int ListAll(ArgList&) const;
    int SetListDebug(ArgList&);
    int ClearList(ArgList&);
    void RemoveDataSet(DataSet*);
    int RemoveDataSet(ArgList&);
    int TrajLength( std::string const&, std::vector<std::string> const&);
    int Run();
    /// Write all DataFiles
    void MasterDataFileWrite();
  private:
    int SetTrajMode(TrajModeType, std::string const&, Topology*, ArgList&);
    int SetTrajMode(TrajModeType);
    int AddReference(DataSet_Coords_REF*, Topology*, DataSet_Coords*, ArgList&, std::string const&, std::string const&, std::string const&);
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

    void ListState() const;
    int RunNormal();
    int RunEnsemble();
#   ifdef MPI
    void DivideFramesAmongProcesses(int&, int&, int&, int, Parallel::Comm const&) const;
    int PreloadCheck(int, int, int&, int&) const;
    int RunParallel();
    int RunParaEnsemble();
    //int RunSingleTrajParallel();
#   endif
    void Init_Timers();
    void Time_Summary() const;
    // -------------------------------------------
    DataSetList DSL_;             ///< List of DataSets
    DataFileList DFL_;            ///< List of DataFiles that DataSets will be written to.
    TrajinList trajinList_;       ///< List of input trajectories/ensembles.
    ActionList actionList_;       ///< List of Actions to be performed during a Run.
    TrajoutList trajoutList_;     ///< List of trajectories to be written during a Run.
    EnsembleOutList ensembleOut_; ///< List of ensembles to be written during a Run.
    AnalysisList analysisList_;   ///< List of Analyses to be performed
    
    int debug_;         ///< Debug level.
    int refDebug_;      ///< Reference debug level.
    int topDebug_;      ///< Topology debug level.
    bool showProgress_; ///< If true, display progress during Run.
    bool quietBlocks_;  ///< If true suppress output when executing control blocks.
    bool exitOnError_;  ///< If true exit when errors encountered instead of continuing.
    bool recordAllInput_; ///< When true save all input to log, even errors.
    /// If true do not process input trajectories when no actions/output trajectories.
    bool noEmptyRun_; // DEBUG: false is used for benchmarking trajectory read speed.
    TrajModeType mode_; ///< Current trajectory mode (NORMAL/ENSEMBLE)
    Timer init_time_;     ///< Run initialization time.
    Timer frames_time_;   ///< Run frame processing time.
    Timer post_time_;     ///< Run post-frame processing (e.g. Action::Print()) time.
    Timer analysis_time_; ///< Run analysis time.
    Timer run_time_;      ///< Total run time.
    Timer write_time_;    ///< Run data file write time.
#   ifdef MPI
    bool forceParallelEnsemble_; ///< If true run parallel ensemble even with 1 thread/member
    Timer sync_time_;     ///< DataSet/Action total sync time.
    Timer master_time_;   ///< Total frame processing time across all ranks.
#   endif
};
#endif
