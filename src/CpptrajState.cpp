#include "CpptrajState.h"
#include "CpptrajStdio.h"
#include "Action_CreateCrd.h" // in case default COORDS need to be created
#include "DataSet_Coords_REF.h" // AddReference
#include "DataSet_Topology.h" // AddTopology
#include "ProgressBar.h"
#ifdef MPI
# include "Parallel.h"
# include "DataSet_Coords_TRJ.h"
# include "EnsembleNavigator.h"
# ifdef TIMER
#   include "EnsembleIn.h"
# endif
#endif

/// CONSTRUCTOR
CpptrajState::CpptrajState() :
  debug_(0),
  refDebug_(0),
  topDebug_(0),
  showProgress_(true),
  quietBlocks_(false),
  exitOnError_(true),
  recordAllInput_(true),
  noEmptyRun_(false),
  mode_(UNDEFINED)
# ifdef MPI
  , forceParallelEnsemble_(false)
# endif
{}

/** This version of SetTrajMode() may be called from AddToActionQueue() to set
  * NORMAL as the default mode without setting up a trajectory. May also be used
  * to clear the mode via UNDEFINED. It is illegal to call this for ENSEMBLE.
  */
int CpptrajState::SetTrajMode(TrajModeType modeIn) {
  if (modeIn == ENSEMBLE) return 1;
  ArgList blankArg;
  return (SetTrajMode( modeIn, std::string(), 0, blankArg ));
}

/** Main routine for setting up trajectory input mode (NORMAL/ENSEMBLE), and
  * potentiall adding a trajectory/ensemble to trajinList_. In parallel this
  * routine is also responsible for initial setup of comms for NORMAL runs.
  */
int CpptrajState::SetTrajMode(TrajModeType modeIn, std::string const& fnameIn,
                              Topology* topIn, ArgList& argIn)
{
  if (modeIn == UNDEFINED) {
    mode_ = modeIn;
    DSL_.SetEnsembleNum( -1 );
    DFL_.SetEnsembleNum( -1 );
#   ifdef MPI
    Parallel::SetupComms( -1 );
#   endif
    return 0;
  }
  if (mode_ != UNDEFINED) {
    if (mode_ != modeIn) {
      mprinterr("Error: 'trajin' and 'ensemble' are mutually exclusive.\n");
      return 1;
    }
  } else // Mode not yet set. Set DataSetList / DataFileList mode if necessary.
    mode_ = modeIn;
  if (mode_ == ENSEMBLE) {
    if (trajinList_.AddEnsembleIn( fnameIn, topIn, argIn )) return 1;
#   ifdef MPI
    // NOTE: SetupComms is called during ensemble setup.
    //rprintf("DEBUG: Inside SetTrajMode(%i): EnsembleComm rank %i\n", (int)modeIn,
    //        Parallel::EnsembleComm().Rank());
    // Make all sets a member of this ensemble.
    DSL_.SetEnsembleNum( Parallel::EnsembleComm().Rank() );
    // This tells all DataFiles to append ensemble member number to file name.
    DFL_.SetEnsembleNum( Parallel::EnsembleComm().Rank() );
#   else
    // Ensemble 0 is set up first. All others are setup in RunEnsemble.
    DSL_.SetEnsembleNum( 0 );
#   endif
  } else if (mode_ == NORMAL) {
    if (topIn != 0) {
      if (trajinList_.AddTrajin( fnameIn, topIn, argIn )) return 1;
    }
#   ifdef MPI
    // NOTE: SetupComms is called here for trajin, whereas in ensemble it is
    //       called inside ensemble setup. This is because ensemble setup
    //       typically requires communication.
    if (Parallel::SetupComms( 1 )) return 1;
#   endif
  }
  return 0;
}

// CpptrajState::AddInputTrajectory()
int CpptrajState::AddInputTrajectory( std::string const& fname ) {
  ArgList args( fname );
  return AddInputTrajectory( args );
}

// CpptrajState::AddInputTrajectory()
int CpptrajState::AddInputTrajectory( ArgList& argIn ) {
  Topology* top = DSL_.GetTopology( argIn );
  if (top == 0) {
    mprinterr("Error: No topology selected or no topologies present.\n");
    return 1;
  }
  return (SetTrajMode( NORMAL, argIn.GetStringNext(), top, argIn ));
}

// CpptrajState::AddInputEnsemble()
int CpptrajState::AddInputEnsemble( ArgList& argIn ) {
  Topology* top = DSL_.GetTopology( argIn );
  if (top == 0) {
    mprinterr("Error: No topology selected or no topologies present.\n");
    return 1;
  }
  return (SetTrajMode( ENSEMBLE, argIn.GetStringNext(), top, argIn ));
}

// CpptrajState::AddOutputTrajectory()
int CpptrajState::AddOutputTrajectory( ArgList& argIn ) {
  // Default to NORMAL if not set.
  if (mode_ == UNDEFINED) {
    mprintf("Warning: Output traj specified before trajin/ensemble. Assuming trajin.\n");
    SetTrajMode( NORMAL );
  }
  std::string fname = argIn.GetStringNext();
  Topology* top = DSL_.GetTopology( argIn );
  int err = 1;
  if (mode_ == NORMAL)
    err = trajoutList_.AddTrajout( fname, argIn, top );
  else if (mode_ == ENSEMBLE)
    err = ensembleOut_.AddEnsembleOut( fname, argIn, top, trajinList_.EnsembleSize() );
  return err;
}

// CpptrajState::AddOutputTrajectory()
int CpptrajState::AddOutputTrajectory( std::string const& fname ) {
  // FIXME Should this use the last Topology instead?
  ArgList tmpArg(fname);
  return AddOutputTrajectory( tmpArg );
}

// CpptrajState::AddToActionQueue()
CpptrajState::RetType CpptrajState::AddToActionQueue( Action* actIn, ArgList& argIn ) {
  argIn.MarkArg(0);
  // Default to NORMAL if not set.
  if (mode_ == UNDEFINED) {
    mprintf("Warning: Action specified before trajin/ensemble. Assuming trajin.\n");
    SetTrajMode( NORMAL );
  }
# ifdef MPI
  DSL_.SetNewSetsNeedSync( true );
  ActionInit init(DSL_, DFL_, Parallel::TrajComm());
# else 
  ActionInit init(DSL_, DFL_);
# endif
  RetType err = OK;
  if (actionList_.AddAction( actIn, argIn, init )) err = ERR;
# ifdef MPI
  DSL_.SetNewSetsNeedSync( false );
  if (Parallel::World().CheckError( err )) err = ERR;
# endif
  return err;
}

// CpptrajState::AddToAnalysisQueue()
CpptrajState::RetType CpptrajState::AddToAnalysisQueue( Analysis* anaIn, ArgList& argIn ) {
  argIn.MarkArg(0);
  AnalysisSetup setup(DSL_, DFL_);
  RetType err = OK;
  if (analysisList_.AddAnalysis( anaIn, argIn, setup )) err = ERR;
# ifdef MPI
  if (Parallel::World().CheckError( err )) err = ERR;
# endif
  return err;
}

// -----------------------------------------------------------------------------
CpptrajState::ListKeyType CpptrajState::ListKeys[] = {
  {L_ACTION,   "actions" }, {L_ACTION,   "action"   },
  {L_TRAJIN,   "trajin"  },
  {L_REF,      "ref"     }, {L_REF,      "reference"},
  {L_TRAJOUT,  "trajout" },
  {L_PARM,     "parm"    }, {L_PARM,     "topology" },
  {L_ANALYSIS, "analysis"}, {L_ANALYSIS, "analyses" },
  {L_DATAFILE, "datafile"}, {L_DATAFILE, "datafiles"},
  {L_DATASET,  "dataset" }, {L_DATASET,  "datasets" }, {L_DATASET, "data"},
  {N_LISTS,    0         }
};

std::string CpptrajState::PrintListKeys() {
  std::string keys;
  for (const ListKeyType* ptr = ListKeys; ptr->Key_ != 0; ptr++) {
    keys += " ";
    keys.append( ptr->Key_ );
  }
  return keys;
}

/** Select lists from ArgList */
std::vector<bool> CpptrajState::ListsFromArg( ArgList& argIn, bool allowEmptyKeyword ) const {
  std::vector<bool> enabled( (int)N_LISTS, false );
  std::string listKeyword = argIn.GetStringNext();
  if (listKeyword.empty() || listKeyword == "all") {
    // See if enabling all lists is permissible.
    if (listKeyword.empty() && !allowEmptyKeyword) {
      mprinterr("Error: A specific list name or 'all' must be specified for '%s'\n",
                argIn.Command());
      return enabled; // All are false
    }
    enabled.assign( (int)N_LISTS, true );
  } else {
    while (!listKeyword.empty()) {
      const ListKeyType* ptr = ListKeys;
      for (; ptr->Key_ != 0; ptr++)
        if (listKeyword == ptr->Key_) {
          enabled[ptr->Type_] = true;
          break;
        }
      if (ptr->Key_ == 0) {
        mprinterr("Error: Unrecognized list name: '%s'\n", listKeyword.c_str());
        mprinterr("Error: Recognized keys:%s\n", PrintListKeys().c_str());
        return std::vector<bool>( (int)N_LISTS, false );
      }
      listKeyword = argIn.GetStringNext();
    }
  }
  return enabled;
}

/** List all members of specified lists */
int CpptrajState::ListAll( ArgList& argIn ) const {
  std::vector<bool> enabled = ListsFromArg( argIn, true );
  if ( enabled[L_ACTION]   ) actionList_.List();
  if ( enabled[L_TRAJIN]   ) trajinList_.List();
  if ( enabled[L_REF]      ) DSL_.ListReferenceFrames();
  if ( enabled[L_TRAJOUT]  ) {
    trajoutList_.List( trajinList_.PindexFrames() );
    ensembleOut_.List( trajinList_.PindexFrames() );
  }
  if ( enabled[L_PARM]     ) DSL_.ListTopologies();
  if ( enabled[L_ANALYSIS] ) analysisList_.List();
  if ( enabled[L_DATAFILE] ) DFL_.List();
  if ( enabled[L_DATASET]  ) DSL_.List();
  return 0;
}

/** Set debug level of specified lists. */
int CpptrajState::SetListDebug( ArgList& argIn ) {
  debug_ = argIn.getNextInteger(0);
  if (debug_ > 0) mprintf("\tGeneral debug level set to %i\n", debug_);
  std::vector<bool> enabled = ListsFromArg( argIn, true );
  if ( enabled[L_ACTION]   ) {
    actionList_.SetDebug( debug_ );
    if (debug_ > 0) mprintf("\tAction debug level set to %i\n", debug_);
  }
  if ( enabled[L_TRAJIN]   ) {
    trajinList_.SetDebug( debug_ );
    if (debug_ > 0) mprintf("\tInput trajectory/ensemble debug level set to %i\n", debug_);
  }
  if ( enabled[L_REF]      ) {
    refDebug_ = debug_;
    if (refDebug_ > 0) mprintf("\tReference debug level set to %i\n", refDebug_);
  }
  if ( enabled[L_TRAJOUT]  ) {
    trajoutList_.SetDebug( debug_ );
    ensembleOut_.SetDebug( debug_ );
    if (debug_ > 0) mprintf("\tOutput trajectory/ensemble debug level set to %i\n", debug_);
  }
  if ( enabled[L_PARM]     ) {
    topDebug_ = debug_;
    if (topDebug_ > 0) mprintf("\tTopology debug level set to %i\n", topDebug_);
  }
  if ( enabled[L_ANALYSIS] ) {
    analysisList_.SetDebug( debug_ );
    if (debug_ > 0) mprintf("\tAnalysis debug level set to %i\n", debug_);
  }
  if ( enabled[L_DATAFILE] ) {
    DFL_.SetDebug( debug_ );
    if (debug_ > 0) mprintf("\tData file debug level set to %i\n", debug_);
  }
  if ( enabled[L_DATASET]  ) {
    DSL_.SetDebug( debug_ );
    if (debug_ > 0) mprintf("\tData set debug level set to %i\n", debug_);
  }
  return 0;
}

/** Clear specified lists */
int CpptrajState::ClearList( ArgList& argIn ) {
  std::vector<bool> enabled = ListsFromArg( argIn, false );
  if ( enabled[L_ACTION]   ) {
    mprintf("\tClearing Actions.\n");
    actionList_.Clear();
  }
  if ( enabled[L_TRAJIN]   ) {
    mprintf("\tClearing input trajectories/ensembles.\n");
    trajinList_.Clear();
    SetTrajMode( UNDEFINED );
  }
  if ( enabled[L_REF]      ) {
    mprintf("\tClearing reference coordinates.\n");
    DSL_.ClearRef();
  }
  if ( enabled[L_TRAJOUT]  ) {
    mprintf("\tClearing output trajectories.\n");
    trajoutList_.Clear();
    ensembleOut_.Clear();
  }
  if ( enabled[L_PARM]     ) {
    mprintf("\tClearing topologies.\n");
    DSL_.ClearTop();
  }
  if ( enabled[L_ANALYSIS] ) {
    mprintf("\tClearing Analyses.\n");
    analysisList_.Clear();
  }
  if ( enabled[L_DATAFILE] ) {
    mprintf("\tClearing data files.\n");
    DFL_.Clear();
  }
  if ( enabled[L_DATASET]  ) {
    mprintf("\tClearing data sets.\n");
    // Make sure sets are cleared from data files as well.
    for (DataSetList::const_iterator ds = DSL_.begin(); ds != DSL_.end(); ++ds)
      DFL_.RemoveDataSet( *ds );
    DSL_.Clear();
  }
  return 0;
}

// CpptrajState::RemoveDataSet()
void CpptrajState::RemoveDataSet(DataSet* dsIn) {
  DFL_.RemoveDataSet( dsIn );
  DSL_.RemoveSet( dsIn );
}

/** Remove DataSet from State */
int CpptrajState::RemoveDataSet( ArgList& argIn ) {
  // Need to first make sure they are removed from DataFiles etc also.
  // FIXME: Currently no good way to check if Actions/Analyses will be
  //        made invalid by DataSet removal.
  std::string removeArg = argIn.GetStringNext();
  if (removeArg.empty()) {
    mprinterr("Error: No data set(s) specified for removal.\n");
    return 1;
  }
  DataSetList tempDSL = DSL_.GetMultipleSets( removeArg );
  if (!tempDSL.empty()) {
    for (DataSetList::const_iterator ds = tempDSL.begin();
                                     ds != tempDSL.end(); ++ds)
    {
      mprintf("\tRemoving \"%s\"\n", (*ds)->legend());
      RemoveDataSet( *ds );
    }
  }
  return 0;
}

// CpptrajState::TrajLength()
// NOTE: MMPBSA.py relies on this.
int CpptrajState::TrajLength( std::string const& topname, 
                              std::vector<std::string> const& trajinFiles)
{
  if (AddTopology( topname, ArgList() )) return 1;
  for (std::vector<std::string>::const_iterator trajinName = trajinFiles.begin();
                                                trajinName != trajinFiles.end();
                                                ++trajinName)
    if (AddInputTrajectory( *trajinName )) return 1;
  loudPrintf("Frames: %i\n", trajinList_.MaxFrames());
  return 0;
}

void CpptrajState::Init_Timers() {
  init_time_.Reset();
  frames_time_.Reset();
  post_time_.Reset();
  analysis_time_.Reset();
  run_time_.Reset();
  write_time_.Reset();
# ifdef MPI
  sync_time_.Reset();
  master_time_.Reset();
# endif
}

void CpptrajState::Time_Summary() const {
  mprintf("\nRUN TIMING:\n");
  init_time_.WriteTiming(2,     "Init               :", run_time_.Total());
# ifdef MPI
  master_time_.WriteTiming(2,   "Trajectory Process :", run_time_.Total());
  sync_time_.WriteTiming(2,     "Data Set Sync      :", run_time_.Total());
# else
  frames_time_.WriteTiming(2,   "Trajectory Process :", run_time_.Total());
# endif
  post_time_.WriteTiming(2,     "Action Post        :", run_time_.Total());
  analysis_time_.WriteTiming(2, "Analysis           :", run_time_.Total());
  write_time_.WriteTiming(2,    "Data File Write    :", run_time_.Total());
  double other_time = run_time_.Total() - init_time_.Total() -
#                     ifdef MPI
                      master_time_.Total() - sync_time_.Total() -
#                     else
                      frames_time_.Total() -
#                     endif
                      post_time_.Total() - analysis_time_.Total() - write_time_.Total();
  mprintf("TIME:\t\tOther              : %.4f s (%6.2f%%)\n",
          other_time, other_time/run_time_.Total());
  run_time_.WriteTiming(1, "Run Total");
}
// -----------------------------------------------------------------------------
// CpptrajState::Run()
int CpptrajState::Run() {
  Init_Timers();
  run_time_.Start();
  int err = 0;
  // Special case: check if _DEFAULTCRD_ COORDS DataSet is defined. If so,
  // this means 1 or more actions has requested that a default COORDS DataSet
  // be created.
  DataSet* default_crd = DSL_.FindSetOfType("_DEFAULTCRD_", DataSet::COORDS);
  if (default_crd != 0) {
    mprintf("Warning: One or more analyses requested creation of default COORDS DataSet.\n");
    // If the DataSet has already been written to do not create again.
    if (default_crd->Size() > 0)
      mprintf("Warning: Default COORDS DataSet has already been written to.\n");
    else {
      // If no input trajectories this will not work.
      if (trajinList_.empty()) {
        mprinterr("Error: Cannot create COORDS DataSet; no input trajectories specified.\n");
        return 1;
      }
#     ifdef MPI
      // Default COORDS DataSet may need to be synced.
      default_crd->SetNeedsSync( true );
#     endif
      ArgList crdcmd("createcrd _DEFAULTCRD_");
      crdcmd.MarkArg(0);
      if (AddToActionQueue( new Action_CreateCrd(), crdcmd ))
        return 1;
    }
  }
  mprintf("---------- RUN BEGIN -------------------------------------------------\n");
  if (trajinList_.empty()) 
    mprintf("Warning: No input trajectories specified.\n");
  else if (actionList_.Empty() && trajoutList_.Empty() && ensembleOut_.Empty() && noEmptyRun_)
    mprintf("Warning: No actions/output trajectories specified.\n");
  else {
    switch ( mode_ ) {
#     ifdef MPI
      case NORMAL:
/*        // TEST - single traj parallel
        if (trajinList_.Size() == 1)
          err = RunSingleTrajParallel();
        else*/
          err = RunParallel();
        break;
      case ENSEMBLE:
        if (Parallel::TrajComm().Size() > 1 || forceParallelEnsemble_)
          err = RunParaEnsemble();
        else
          err = RunEnsemble();
        break;
#     else
      case NORMAL   : err = RunNormal(); break;
      case ENSEMBLE : err = RunEnsemble(); break;
#     endif
      case UNDEFINED: break;
    }
    // Clean up Actions if run completed successfully.
    if (err == 0) {
      actionList_.Clear();
      trajoutList_.Clear();
      ensembleOut_.Clear();
      DSL_.SetDataSetsPending(false);
    }
  }
  // Run Analyses if any are specified.
  if (err == 0)
    err = RunAnalyses();
  write_time_.Start();
  if (err == 0 || !exitOnError_) {
    DSL_.ListDataOnly();
    // Print DataFile information and write DataFiles
    DFL_.List();
    MasterDataFileWrite();
  }
  write_time_.Stop();
  run_time_.Stop();
  Time_Summary();
  mprintf("---------- RUN END ---------------------------------------------------\n");
  return err;
}

/** Prior to a run, list the current state. */
void CpptrajState::ListState() const {
  if (mode_ == ENSEMBLE) trajinList_.FirstEnsembleReplicaInfo();
  DSL_.ListTopologies();
  trajinList_.List();
  DSL_.ListReferenceFrames();
  if (mode_ == ENSEMBLE)
    ensembleOut_.List( trajinList_.PindexFrames() );
  else
    trajoutList_.List( trajinList_.PindexFrames() );
}

// -----------------------------------------------------------------------------
// CpptrajState::RunEnsemble()
/** Process ensemble in serial or parallel, individual trajectories in serial. */
int CpptrajState::RunEnsemble() {
  init_time_.Start();
  // Print ensemble size
  int ensembleSize = trajinList_.EnsembleSize();
  mprintf("\nENSEMBLE INFO:\n  Ensemble size is %i\n", ensembleSize);
  // Allocate space to hold position of each incoming frame in replica space.
# ifdef MPI
  // In parallel only two frames needed; one for reading, one for receiving.
  FramePtrArray SortedFrames( 2 );
  FrameArray FrameEnsemble( 2 );
  // Each thread will process one member of the ensemble, so local ensemble
  // size is effectively 1.
  ensembleSize = 1;
# else
  FramePtrArray SortedFrames( ensembleSize );
  FrameArray FrameEnsemble( ensembleSize );
# endif
  // List state
  ListState();

  // Allocate DataSets in the master DataSetList based on # frames to be read
  DSL_.AllocateSets( trajinList_.MaxFrames() );

  // Allocate an ActionList for each member of the ensemble.
  std::vector<ActionList*> ActionEnsemble( ensembleSize );
  ActionEnsemble[0] = &actionList_;
  for (int member = 1; member < ensembleSize; member++)
    ActionEnsemble[member] = new ActionList();
  // If we are on a single thread, give each member its own copy of the
  // current topology address. This way if topology is modified by a member,
  // e.g. in strip or closest, subsequent members wont be trying to modify 
  // an already-modified topology.
  std::vector<ActionSetup> EnsembleParm( ensembleSize );
  // Give each member its own copy of current frame address. This way if 
  // frame is modified by a member things like trajout know about it.
  FramePtrArray CurrentFrames( ensembleSize );
# ifndef MPI
  // Silence action output for members > 0.
  if (debug_ == 0) SetWorldSilent( true ); 
  // Set up Actions for each ensemble member > 0.
  for (int member = 1; member < ensembleSize; ++member) {
    // All DataSets that will be set up will be part of this ensemble 
    DSL_.SetEnsembleNum( member );
    ActionInit init( DSL_, DFL_ );
    // Initialize actions for this ensemble member based on original actionList_
    if (!actionList_.Empty()) {
      if (debug_ > 0) mprintf("***** ACTIONS FOR ENSEMBLE MEMBER %i:\n", member);
      for (int iaction = 0; iaction < actionList_.Naction(); iaction++) { 
        ArgList actionArgs = actionList_.ActionArgs(iaction);
        // Attempt to add same action to this ensemble. 
        if (ActionEnsemble[member]->AddAction(actionList_.ActionAlloc(iaction), actionArgs, init))
            return 1;
      }
    }
  }
  SetWorldSilent( false );
# endif
  init_time_.Stop();

  // ========== A C T I O N  P H A S E ==========
  int actionSet = 0;            // Internal data frame
  int readSets = 0;             // Number of frames actually read
  int lastPindex = -1;          // Index of the last loaded parm file
  CoordinateInfo lastCoordInfo; // Coordinate info for previous ensemble
  // Loop over every ensemble 
# ifdef TIMER
  Timer trajin_time;
  Timer setup_time;
  Timer actions_time;
  Timer trajout_time;
# endif
  frames_time_.Start();
  mprintf("\nBEGIN ENSEMBLE PROCESSING:\n");
  ProgressBar progress;
  if (showProgress_)
    progress.SetupProgress( trajinList_.MaxFrames() );
  for ( TrajinList::ensemble_it ens = trajinList_.ensemble_begin();
                                ens != trajinList_.ensemble_end(); ++ens)
  {
    // Open up the ensemble. If an error occurs, bail
    if ( (*ens)->BeginEnsemble() ) {
      mprinterr("Error: Could not open ensemble %s.\n",(*ens)->Traj().Filename().full());
      break;
    }
    // Set current parm from current ensemble.
    Topology* currentTop = (*ens)->Traj().Parm();
    CoordinateInfo const& currentCoordInfo = (*ens)->EnsembleCoordInfo();
    currentTop->SetBoxFromTraj( currentCoordInfo.TrajBox() ); // FIXME necessary?
    int topFrames = trajinList_.TopFrames( currentTop->Pindex() );
    for (int member = 0; member < ensembleSize; ++member)
      EnsembleParm[member].Set( currentTop, currentCoordInfo, topFrames );
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != currentTop->Pindex());
#   ifdef TIMER
    setup_time.Start();
#   endif
    // If Parm has changed or trajectory frame has changed, reset the frame. 
    if (parmHasChanged || currentCoordInfo != lastCoordInfo) {
      FrameEnsemble.SetupFrames(currentTop->Atoms(), currentCoordInfo);
      lastCoordInfo = currentCoordInfo;
    }
    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set up actions for this parm
      bool setupOK = true;
      for (int member = 0; member < ensembleSize; ++member) {
        // Silence action output for all beyond first member.
        if (member > 0)
          SetWorldSilent( true );
        if (ActionEnsemble[member]->SetupActions( EnsembleParm[member], exitOnError_ )) {
#         ifdef MPI
          rprintf("Warning: Ensemble member %i: Could not set up actions for %s: skipping.\n",
                  Parallel::EnsembleComm().Rank(), EnsembleParm[member].Top().c_str());
#         else
          mprintf("Warning: Ensemble member %i: Could not set up actions for %s: skipping.\n",
                  member, EnsembleParm[member].Top().c_str());
#         endif
          setupOK = false;
        }
      }
      // Re-enable output
      SetWorldSilent( false );
      if (!setupOK) continue;
      // Set up any related output ensembles.
      // TODO: Currently assuming topology is always modified the same
      //       way for all actions. If this behavior ever changes the
      //       following line will cause undesireable behavior.
      ensembleOut_.SetupEnsembleOut( EnsembleParm[0].TopAddress(),
                                     EnsembleParm[0].CoordInfo(),
                                     EnsembleParm[0].Nframes() );
      lastPindex = currentTop->Pindex();
    }
#   ifdef TIMER
    setup_time.Stop();
#   endif
    // Loop over every collection of frames in the ensemble
    (*ens)->Traj().PrintInfoLine();
#   ifdef TIMER
    trajin_time.Start();
    bool readMoreFrames = (*ens)->GetNextEnsemble(FrameEnsemble, SortedFrames);
    trajin_time.Stop();
    while ( readMoreFrames )
#   else
    while ( (*ens)->GetNextEnsemble(FrameEnsemble, SortedFrames) )
#   endif
    {
      if (!(*ens)->BadEnsemble()) {
        bool suppress_output = false;
        for (int member = 0; member != ensembleSize; ++member) {
          // Since Frame can be modified by actions, save original and use currentFrame
          ActionFrame currentFrame( SortedFrames[member], actionSet );
          //rprintf("DEBUG: currentFrame=%x SortedFrames[0]=%x\n",currentFrame, SortedFrames[0]);
          if ( currentFrame.Frm().CheckCoordsInvalid() )
            rprintf("Warning: Ensemble member %i frame %i may be corrupt.\n",
                    member, (*ens)->Traj().Counter().PreviousFrameNumber()+1);
#         ifdef TIMER
          actions_time.Start();
#         endif
          // Perform Actions on Frame
          suppress_output = ActionEnsemble[member]->DoActions(actionSet, currentFrame);
          CurrentFrames[member] = currentFrame.FramePtr();
#         ifdef TIMER
          actions_time.Stop();
#         endif
        } // END loop over actions
        // Do Output
        if (!suppress_output) {
#         ifdef TIMER
          trajout_time.Start();
#         endif 
          if (ensembleOut_.WriteEnsembleOut(actionSet, CurrentFrames))
          {
            mprinterr("Error: Writing ensemble output traj, frame %i\n", actionSet+1);
            if (exitOnError_) return 1; 
          }
#         ifdef TIMER
          trajout_time.Stop();
#         endif
        }
      } else {
#       ifdef MPI
        rprinterr("Error: Could not read frame %i for ensemble.\n", actionSet + 1);
#       else
        mprinterr("Error: Could not read frame %i for ensemble.\n", actionSet + 1);
#       endif
      }
      if (showProgress_) progress.Update( actionSet );
      // Increment frame counter
      ++actionSet;
#     ifdef TIMER
      trajin_time.Start();
      readMoreFrames = (*ens)->GetNextEnsemble(FrameEnsemble, SortedFrames);
      trajin_time.Stop();
#     endif
    }

    // Close the trajectory file
    (*ens)->EndEnsemble();
    // Update how many frames have been processed.
    readSets += (*ens)->Traj().Counter().NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  mprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);
  frames_time_.Stop();
  mprintf("TIME: Avg. throughput= %.4f frames / second.\n",
          (double)readSets / frames_time_.Total());
# ifdef MPI
  master_time_ = frames_time_;
# endif
# ifdef TIMER
  trajin_time.WriteTiming(1,  "Trajectory read:        ", frames_time_.Total());
  setup_time.WriteTiming(1,   "Action setup:           ", frames_time_.Total());
  actions_time.WriteTiming(1, "Action frame processing:", frames_time_.Total());
  trajout_time.WriteTiming(1, "Trajectory output:      ", frames_time_.Total());
# ifdef MPI
  EnsembleIn::TimingData(trajin_time.Total());
# endif
# endif
  // Close output trajectories
  ensembleOut_.CloseEnsembleOut();
  post_time_.Start();
  // ========== A C T I O N  O U T P U T  P H A S E ==========
  mprintf("\nENSEMBLE ACTION OUTPUT:\n");
  for (int member = 0; member < ensembleSize; ++member)
    ActionEnsemble[member]->PrintActions();
  post_time_.Stop();
  // Clean up ensemble action lists
  for (int member = 1; member < ensembleSize; member++)
    delete ActionEnsemble[member];

  return 0;
}
#ifdef MPI
// -----------------------------------------------------------------------------
void CpptrajState::DivideFramesAmongThreads(int& my_start, int& my_stop, int& my_frames,
                                            int maxFrames, Parallel::Comm const& commIn) const
{
  my_frames = commIn.DivideAmongThreads(my_start, my_stop, maxFrames);
  std::vector<int> frames_per_thread( commIn.Size() );
  commIn.GatherMaster(&my_frames, 1, MPI_INT, &frames_per_thread[0]);
  // Print how many frames each rank will process.
  if (commIn.Master()) {
    mprintf("\nPARALLEL INFO:\n");
    if (Parallel::EnsembleComm().Size() > 1)
      mprintf("  %i threads per ensemble member.\n", commIn.Size());
    for (int rank = 0; rank != commIn.Size(); rank++)
      mprintf("  Thread %i will process %i frames.\n", rank, frames_per_thread[rank]);
  }
  commIn.Barrier();
  if (debug_ > 0) rprintf("Start %i Stop %i Frames %i\n", my_start+1, my_stop, my_frames);
}

/** Figure out if any frames need to be preloaded on ranks. Should NOT be 
  * called by master.
  */
int CpptrajState::PreloadCheck(int my_start, int my_frames,
                               int& n_previous_frames, int& preload_start) const
{
  n_previous_frames = actionList_.NumPreviousFramesReqd();
  if (n_previous_frames < 1) return 0;
  preload_start = my_start - n_previous_frames;
  if (debug_ > 0)
    rprintf("DEBUG: Preloading frames from %i to %i\n", my_start - n_previous_frames, my_start-1);
  if (preload_start < 0) {
    rprinterr("Error: Cannot preload, start is before beginning of traj.\n");
    return 1;
  } else if (my_frames == n_previous_frames) {
    rprinterr("Error: Number of preload frames is same as number of processed frames.\n"
              "Error:   Try reducing the number of threads.\n");
    return 1;
  } 
  if (n_previous_frames > (my_frames / 2))
    rprintf("Warning: Number of preload frames is greater than half the "
                      "number of processed frames.\n"
            "Warning:   Try reducing the number of threads.\n");
  rprintf("Warning: Preloading %i frames. These frames will NOT have Actions performed on them.\n",
          n_previous_frames);
  return 0;
}
// -----------------------------------------------------------------------------
/** Process ensemble and individual trajectories in parallel. */
int CpptrajState::RunParaEnsemble() {
  init_time_.Start();
  // Set Comms
  Parallel::Comm const& EnsComm = Parallel::EnsembleComm();
  Parallel::Comm const& TrajComm = Parallel::TrajComm();
  // Print ensemble size
  mprintf("\nENSEMBLE INFO:\n  Ensemble size is %i\n", trajinList_.EnsembleSize());
  // Only two frames needed; one for reading, one for receiving.
  FrameArray FrameEnsemble( 2 );
  FramePtrArray SortedFrames( 2 );
  // In parallel, each ensemble member is responsible for only 1 frame.
  FramePtrArray CurrentFrames( 1 );
  // List state
  ListState();
  // Currently only let this happen if all trajectories share same topology.
  EnsembleNavigator NAV;
  int err = NAV.AddEnsembles(trajinList_.ensemble_begin(), trajinList_.ensemble_end());
  if (Parallel::World().CheckError( err )) return 1;

  // Divide frames among threads
  int my_start, my_stop, my_frames;
  DivideFramesAmongThreads(my_start, my_stop, my_frames, NAV.IDX().MaxFrames(), TrajComm);
  // Ensure at least 1 frame per thread, otherwise some ranks could cause hangups.
  if (my_frames > 0)
    err = 0;
  else {
    rprinterr("Error: Thread is processing less than 1 frame. Try reducing # threads.\n");
    err = 1;
  }
  if (Parallel::World().CheckError( err )) return 1;

  // Allocate DataSets in the master DataSetList based on # frames to be read by this thread.
  DSL_.AllocateSets( my_frames );
  // Any DataSets added to the DataSetList during run will need to be synced.
  DSL_.SetNewSetsNeedSync( true );

  // ----- SETUP PHASE ---------------------------
  NAV.FirstParm()->SetBoxFromTraj( NAV.EnsCoordInfo().TrajBox() ); // FIXME necessary?
  ActionSetup currentParm( NAV.FirstParm(), NAV.EnsCoordInfo(), NAV.IDX().MaxFrames() );
  err = actionList_.SetupActions( currentParm, exitOnError_ );
  if (Parallel::World().CheckError( err )) {
    mprinterr("Error: Could not set up actions for '%s'\n", NAV.FirstParm()->c_str());
    return 1;
  }
  // Set up any related output trajectories.
  if (ensembleOut_.ParallelSetupEnsembleOut( currentParm.TopAddress(),
                                             currentParm.CoordInfo(),
                                             currentParm.Nframes(), TrajComm ))
    return 1;
  // Allocate FrameEnsemble here in case preload is needed.
  FrameEnsemble.SetupFrames( NAV.FirstParm()->Atoms(), NAV.EnsCoordInfo() );
  // Figure out if any frames need to be preloaded on ranks
  int preload_err = 0;
  if (!TrajComm.Master()) {
    int n_previous_frames, preload_start;
    preload_err = PreloadCheck( my_start, my_frames, n_previous_frames, preload_start );
    if (n_previous_frames > 0 && preload_err == 0) {
      Action::FArray preload_frames;
      preload_frames.reserve( n_previous_frames );
      int idx = 0;
      for (int set = preload_start; set != my_start; set++, idx++) {
        NAV.GetEnsemble( set, FrameEnsemble, SortedFrames );
        preload_frames.push_back( *(SortedFrames[0]) );
      }
      preload_err = actionList_.ParallelProcessPreload( preload_frames );
    }
  }
  if (Parallel::World().CheckError( preload_err ) && exitOnError_) return 1;
  init_time_.Stop();
  // ----- ACTION PHASE --------------------------
  frames_time_.Start();
  master_time_.Start();
  mprintf("\nBEGIN PARALLEL ENSEMBLE PROCESSING:\n");
  ProgressBar progress;
  if (showProgress_)
    progress.SetupProgress( my_frames );
  int actionSet = 0; // Internal data frame
  err = 0;
  for (int set = my_start; set != my_stop; set++, actionSet++) {
    // Read the ensemble
    err = NAV.GetEnsemble( set, FrameEnsemble, SortedFrames ); // TODO better error check
    if (err != 0) { 
      rprinterr("Error: Could not open ensemble %i '%s'\n", NAV.IDX().CurrentTrajNum(),
                NAV.CurrentEns()->Traj().Filename().full());
      break;
    }
    if (!NAV.CurrentEns()->BadEnsemble()) {
      // Since Frame can be modified by actions, save original and use currentFrame
      ActionFrame currentFrame( SortedFrames[0], set );
      if ( currentFrame.Frm().CheckCoordsInvalid() )
        rprintf("Warning: Ensemble member %i frame %i may be corrupt.\n",
                EnsComm.Rank(), NAV.CurrentEns()->Traj().Counter().PreviousFrameNumber()+1);
      // Perform Actions on Frame
      bool suppress_output = actionList_.DoActions(actionSet, currentFrame);
      if (!suppress_output) {
        CurrentFrames[0] = currentFrame.FramePtr();
        if (ensembleOut_.WriteEnsembleOut(set, CurrentFrames))
        {
          rprinterr("Error: Writing ensemble output traj, frame %i\n", actionSet+1);
          if (exitOnError_) return 1;
        }
      }
    } else {
      rprinterr("Error: Could not read set %i for ensemble.\n", set + 1);
    }
    if (showProgress_) progress.Update( actionSet );
  } // END loop over frames
  frames_time_.Stop();
  // Close the trajectory file
  NAV.CurrentEns()->EndEnsemble();
  // Collect FPS stats from each rank. 
  std::vector<double> darray( Parallel::World().Size() );
  double rank_fps = (double)actionSet / frames_time_.Total();
  Parallel::World().GatherMaster( &rank_fps, 1, MPI_DOUBLE, &darray[0] ); // Acts as barrier
  master_time_.Stop();
  int global = 0;
  for (int member = 0; member < EnsComm.Size(); member++) {
    mprintf("TIME: Member %4i ranks FPS=", member);
    for (int rank = 0; rank < TrajComm.Size(); rank++)
      mprintf(" %8.2f", darray[global++]);
    mprintf("\n");
  }
  mprintf("TIME: Avg. throughput= %.4f frames / second.\n",
          (double)NAV.IDX().MaxFrames() / master_time_.Total());
# ifdef TIMER
  EnsembleIn::TimingData(master_time_.Total());
# endif
  // Close output trajectories
  ensembleOut_.CloseEnsembleOut();
  DSL_.SetNewSetsNeedSync( false );
  sync_time_.Start();
  // Sync Actions to master thread
  actionList_.SyncActions();
  // Sync data sets to master thread
  if (DSL_.SynchronizeData( TrajComm )) return 1;
  sync_time_.Stop();
  mprintf("\nACTION OUTPUT:\n");
  post_time_.Start();
  // Only call print for master
  if (TrajComm.Master())
    actionList_.PrintActions();
  TrajComm.Barrier();
  post_time_.Stop();
  return 0;
}

// -----------------------------------------------------------------------------
/** Process trajectories in trajinList in parallel. */
int CpptrajState::RunParallel() {
  init_time_.Start();
  // Set comms
  Parallel::Comm const& TrajComm = Parallel::TrajComm();
  // List state 
  ListState();

  // Currently only let this happen if all trajectories share same topology.
  Topology* FirstParm = 0;
  int err = 0;
  for ( TrajinList::trajin_it traj = trajinList_.trajin_begin();
                              traj != trajinList_.trajin_end(); ++traj)
    if (FirstParm == 0)
      FirstParm = (*traj)->Traj().Parm();
    else {
      if (FirstParm != (*traj)->Traj().Parm()) {
        err = 1;
        break;
      }
  }
  if (TrajComm.CheckError( err )) {
    mprinterr("Error: 'trajin' in parallel currently requires all trajectories use"
              " the same topology file.\n");
    return 1;
  }

  // Put all trajectories into a DataSet_Coords_TRJ for random access.
  // FIXME error check above goes here instead?
  DataSet_Coords_TRJ input_traj;
  for ( TrajinList::trajin_it traj = trajinList_.trajin_begin();
                              traj != trajinList_.trajin_end(); ++traj)
    if (input_traj.AddInputTraj( *traj )) { err = 1; break; }
  if (TrajComm.CheckError( err )) return 1;

  // Divide frames among threads.
  int my_start, my_stop, my_frames;
  DivideFramesAmongThreads(my_start, my_stop, my_frames, input_traj.Size(), TrajComm);
  // Ensure at least 1 frame per thread, otherwise some ranks could cause hangups.
  if (my_frames > 0)
    err = 0;
  else {
    rprinterr("Error: Thread is processing less than 1 frame. Try reducing # threads.\n");
    err = 1;
  }
  if (TrajComm.CheckError( err )) return 1;

  // Allocate DataSets in DataSetList based on # frames read by this thread.
  DSL_.AllocateSets( my_frames );
  // Any DataSets added to the DataSetList during run will need to be synced.
  DSL_.SetNewSetsNeedSync( true );

  // ----- SETUP PHASE ---------------------------
  CoordinateInfo const& currentCoordInfo = input_traj.CoordsInfo();
  Topology* top = input_traj.TopPtr();
  top->SetBoxFromTraj( currentCoordInfo.TrajBox() ); // FIXME necessary?
  int topFrames = trajinList_.TopFrames( top->Pindex() );
  ActionSetup currentParm( top, currentCoordInfo, topFrames );
  err = actionList_.SetupActions( currentParm, exitOnError_ );
  if (TrajComm.CheckError( err )) {
    mprinterr("Error: Could not set up actions for '%s'\n", top->c_str());
    return 1;
  }

  // Set up any related output trajectories. 
  if (trajoutList_.ParallelSetupTrajout( currentParm.TopAddress(), currentParm.CoordInfo(),
                                         input_traj.Size(), TrajComm ))
    return 1;

  // Figure out if any frames need to be preloaded on ranks
  int preload_err = 0;
  if (!TrajComm.Master()) {
    int n_previous_frames, preload_start;
    preload_err = PreloadCheck( my_start, my_frames, n_previous_frames, preload_start );
    if (n_previous_frames > 0 && preload_err == 0) {
      Action::FArray preload_frames( n_previous_frames, input_traj.AllocateFrame() );
      int idx = 0;
      for (int set = preload_start; set != my_start; set++, idx++)
        input_traj.GetFrame(set, preload_frames[idx]);
      preload_err = actionList_.ParallelProcessPreload( preload_frames );
    }
  }
  if (TrajComm.CheckError( preload_err ) && exitOnError_) return 1;
  init_time_.Stop();
  // ----- ACTION PHASE --------------------------
  mprintf("\nBEGIN PARALLEL TRAJECTORY PROCESSING:\n");
  frames_time_.Start();
  master_time_.Start();
  ProgressBar progress;
  if (showProgress_)
    progress.SetupProgress( my_frames );
  Frame TrajFrame = input_traj.AllocateFrame();
  int actionSet = 0; // Internal data frame
  for (int set = my_start; set != my_stop; set++, actionSet++) {
    input_traj.GetFrame(set, TrajFrame);
    if (TrajFrame.CheckCoordsInvalid()) // TODO actual frame #
      rprintf("Warning: Set %i coords 1 & 2 overlap at origin; may be corrupt.\n", set + 1);
    ActionFrame currentFrame( &TrajFrame, set );
    bool suppress_output = actionList_.DoActions(actionSet, currentFrame);
    // Trajectory output
    if (!suppress_output) {
      if (trajoutList_.WriteTrajout(set, currentFrame.Frm())) {
        if (exitOnError_) return 1;
      }
    }
    if (showProgress_) progress.Update( actionSet );
  }
  frames_time_.Stop();
  // Collect FPS stats from each rank.
  std::vector<double> darray( TrajComm.Size() );
  double rank_fps = (double)actionSet / frames_time_.Total();
  TrajComm.GatherMaster( &rank_fps, 1, MPI_DOUBLE, &darray[0] ); // Acts as barrier
  master_time_.Stop();
  for (int rank = 0; rank < TrajComm.Size(); rank++)
    mprintf("TIME: Rank %i throughput= %.4f frames / second.\n", rank, darray[rank]);
  mprintf("TIME: Avg. throughput= %.4f frames / second.\n",
          (double)input_traj.Size() / master_time_.Total());
  trajoutList_.CloseTrajout();
  DSL_.SetNewSetsNeedSync( false );
  sync_time_.Start();
  // Sync Actions to master thread
  actionList_.SyncActions();
  // Sync data sets to master thread
  if (DSL_.SynchronizeData( TrajComm )) return 1;
  sync_time_.Stop();
  post_time_.Start();
  mprintf("\nACTION OUTPUT:\n");
  // Only call print for master
  if (TrajComm.Master())
    actionList_.PrintActions();
  TrajComm.Barrier();
  post_time_.Stop();
  return 0;
}

// -----------------------------------------------------------------------------
/*
int CpptrajState::RunSingleTrajParallel() {
  mprintf("DEBUG: Experimental: Opening single trajectory in parallel.\n");
  // Print information.
  DSL_.ListTopologies();
  trajinList_.List();
  DSL_.ListReferenceFrames();
  trajoutList_.List( trajinList_.PindexFrames() );
  // Set up single trajectory for parallel read.
  Trajin* trajin = *(trajinList_.trajin_begin());
  trajin->ParallelBeginTraj( Parallel::World() );
  // Divide frames among threads.
  int total_read_frames = trajin->Traj().Counter().TotalReadFrames();
  int my_start, my_stop, my_frames;
  std::vector<int> rank_frames = DivideFramesAmongThreads(my_start, my_stop, my_frames,
                                                          total_read_frames,
                                                          Parallel::World().Size(),
                                                          Parallel::World().Rank(),
                                                          Parallel::World().Master());
  // Update my start and stop based on offset.
  int traj_offset = trajin->Traj().Counter().Offset();
  int traj_start  = my_start * traj_offset + trajin->Traj().Counter().Start();
  int traj_stop   = my_stop  * traj_offset + trajin->Traj().Counter().Start();
  rprintf("Start and stop adjusted for offset: %i to %i\n", traj_start, traj_stop);
  Parallel::World().Barrier();

  // Allocate DataSets in DataSetList based on # frames read by this thread.
  DSL_.AllocateSets( my_frames );

  // ----- SETUP PHASE ---------------------------
  CoordinateInfo const& currentCoordInfo = trajin->TrajCoordInfo();
  Topology* top = trajin->Traj().Parm();
  top->SetBoxFromTraj( currentCoordInfo.TrajBox() ); // FIXME necessary?
  int topFrames = trajinList_.TopFrames( top->Pindex() );
  ActionSetup currentParm( top, currentCoordInfo, topFrames );
  int err = actionList_.SetupActions( currentParm, exitOnError_ );
  if (Parallel::World().CheckError( err )) {
    mprinterr("Error: Could not set up actions for '%s'\n", top->c_str());
    return 1;
  }

  // Set up any related output trajectories. 
  if (trajoutList_.ParallelSetupTrajout( top, currentCoordInfo,
                                         total_read_frames, Parallel::World() ))
    return 1;

  // ----- ACTION PHASE --------------------------
  Timer frames_time;
  frames_time.Start();
  ProgressBar progress;
  if (showProgress_)
    progress.SetupProgress( my_frames );
  Frame TrajFrame;
  TrajFrame.SetupFrameV(top->Atoms(), currentCoordInfo);
  int actionSet = 0; // Internal data frame
  int trajoutSet = my_start; // Trajout index
  for (int trajinSet = traj_start; trajinSet < traj_stop; trajinSet += traj_offset,
                                                          actionSet++, trajoutSet++)
  {
    trajin->ParallelReadTrajFrame(trajinSet, TrajFrame);
    if (TrajFrame.CheckCoordsInvalid())
      rprintf("Warning: Frame %i coords 1 & 2 overlap at origin; may be corrupt.\n", trajinSet + 1);
    ActionFrame currentFrame( &TrajFrame );
    bool suppress_output = actionList_.DoActions(actionSet, currentFrame);
    // Trajectory output
    if (!suppress_output) {
      if (trajoutList_.WriteTrajout(trajoutSet, currentFrame.Frm())) {
        if (exitOnError_) return 1;
      }
    }
    if (showProgress_) progress.Update( actionSet );
  }
  frames_time.Stop();
  rprintf("TIME: Avg. throughput= %.4f frames / second.\n",
          (double)actionSet / frames_time.Total());
  Parallel::World().Barrier();
  trajin->ParallelEndTraj();
  mprintf("TIME: Avg. throughput= %.4f frames / second.\n",
          (double)total_read_frames / frames_time.Total());
  trajoutList_.CloseTrajout();
  // Sync data sets to master thread
  Timer time_sync;
  time_sync.Start();
  if (DSL_.SynchronizeData( total_read_frames, rank_frames, Parallel::World() )) return 1;
  // Sync Actions to master thread
  actionList_.SyncActions();
  time_sync.Stop();
  time_sync.WriteTiming(1, "Data set/actions sync");
  mprintf("\nACTION OUTPUT:\n");
  // Only call print for master
  if (Parallel::World().Master())
    actionList_.PrintActions();
  Parallel::World().Barrier();
  return 0;
}
*/
#endif
// -----------------------------------------------------------------------------
// CpptrajState::RunNormal()
/** Process trajectories in trajinList. Each frame in trajinList is sent
 *  to the actions in actionList for processing.
 */
int CpptrajState::RunNormal() {
  init_time_.Start();
  // List state
  ListState();
  // Allocate DataSets in the master DataSetList based on # frames to be read
  DSL_.AllocateSets( trajinList_.MaxFrames() );
  init_time_.Stop();
  
  // ========== A C T I O N  P H A S E ==========
  int actionSet = 0;            // Internal data frame
  int readSets = 0;             // Number of frames actually read
  int lastPindex = -1;          // Index of the last loaded parm file
  CoordinateInfo lastCoordInfo; // Coordinate info for previous trajectory.
  Frame TrajFrame;              // Original Frame read in from traj
  // Loop over every trajectory in trajFileList
# ifdef TIMER
  Timer trajin_time;
  Timer setup_time;
  Timer actions_time;
  Timer trajout_time;
# endif
  frames_time_.Start();
  mprintf("\nBEGIN TRAJECTORY PROCESSING:\n");
  ProgressBar progress;
  if (showProgress_)
    progress.SetupProgress( trajinList_.MaxFrames() );
  for ( TrajinList::trajin_it traj = trajinList_.trajin_begin();
                              traj != trajinList_.trajin_end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail 
    if ( (*traj)->BeginTraj() ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->Traj().Filename().full());
      break;
    }
    // Set current parm from current traj.
    Topology* top = (*traj)->Traj().Parm();
    top->SetBoxFromTraj( (*traj)->TrajCoordInfo().TrajBox() ); // FIXME necessary?
    ActionSetup currentSetup( top, (*traj)->TrajCoordInfo(),
                             trajinList_.TopFrames( top->Pindex() ) );
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != currentSetup.Top().Pindex());
#   ifdef TIMER
    setup_time.Start();
#   endif
    // If Parm has changed or trajectory frame has changed, reset the frame.
    if (parmHasChanged || currentSetup.CoordInfo() != lastCoordInfo) {
      TrajFrame.SetupFrameV(currentSetup.Top().Atoms(), currentSetup.CoordInfo());
      lastCoordInfo = currentSetup.CoordInfo();
    }
    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set up actions for this parm
      if (actionList_.SetupActions( currentSetup, exitOnError_ )) {
        mprintf("Warning: Could not set up actions for %s: skipping.\n",
                currentSetup.Top().c_str());
        continue;
      }
      // Set up any related output trajectories 
      trajoutList_.SetupTrajout( currentSetup.TopAddress(),
                                 currentSetup.CoordInfo(),
                                 currentSetup.Nframes() );
      lastPindex = currentSetup.Top().Pindex();
    }
#   ifdef TIMER
    setup_time.Stop();
#   endif
    // Loop over every Frame in trajectory
    (*traj)->Traj().PrintInfoLine();
#   ifdef TIMER
    trajin_time.Start();
    bool readMoreFrames = (*traj)->GetNextFrame(TrajFrame);
    trajin_time.Stop();
    while ( readMoreFrames )
#   else
    while ( (*traj)->GetNextFrame(TrajFrame) )
#   endif
    {
      // Since Frame can be modified by actions, save original and use currentFrame
      ActionFrame currentFrame( &TrajFrame, actionSet );
      // Check that coords are valid.
      if ( currentFrame.Frm().CheckCoordsInvalid() )
        mprintf("Warning: Frame %i coords 1 & 2 overlap at origin; may be corrupt.\n",
                (*traj)->Traj().Counter().PreviousFrameNumber()+1);
        // Perform Actions on Frame
#       ifdef TIMER
        actions_time.Start();
#       endif
        bool suppress_output = actionList_.DoActions(actionSet, currentFrame);
#       ifdef TIMER
        actions_time.Stop();
#       endif
        // Do Output
        if (!suppress_output) {
#         ifdef TIMER
          trajout_time.Start();
#         endif
          if (trajoutList_.WriteTrajout(actionSet, currentFrame.Frm())) {
            if (exitOnError_) return 1;
          }
#         ifdef TIMER
          trajout_time.Stop();
#         endif
        }
      if (showProgress_) progress.Update( actionSet );
      // Increment frame counter
      ++actionSet;
#     ifdef TIMER
      trajin_time.Start();
      readMoreFrames = (*traj)->GetNextFrame(TrajFrame);
      trajin_time.Stop();
#     endif 
    }

    // Close the trajectory file
    (*traj)->EndTraj();
    // Update how many frames have been processed.
    readSets += (*traj)->Traj().Counter().NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  mprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);
  frames_time_.Stop();
  mprintf("TIME: Avg. throughput= %.4f frames / second.\n", 
          (double)readSets / frames_time_.Total());
# ifdef TIMER
  DSL_.Timing();
  trajin_time.WriteTiming(1,  "Trajectory read:        ", frames_time_.Total());
  setup_time.WriteTiming(1,   "Action setup:           ", frames_time_.Total());
  actions_time.WriteTiming(1, "Action frame processing:", frames_time_.Total());
  trajout_time.WriteTiming(1, "Trajectory output:      ", frames_time_.Total());
# endif
  // Close output trajectories.
  trajoutList_.CloseTrajout();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  mprintf("\nACTION OUTPUT:\n");
  post_time_.Start();
  actionList_.PrintActions();
  post_time_.Stop();
  return 0;
}

// -----------------------------------------------------------------------------
// CpptrajState::MasterDataFileWrite()
/** Trigger write of all pending DataFiles. When in parallel ensemble mode,
  * each member of the ensemble will write data to separate files with 
  * numeric extensions.
  */
void CpptrajState::MasterDataFileWrite() { DFL_.WriteAllDF(); }

// CpptrajState::RunAnalyses()
int CpptrajState::RunAnalyses() {
  if (analysisList_.Empty()) return 0;
  analysis_time_.Reset();
  analysis_time_.Start();
  int err = 0;
# ifdef MPI
  // Only master performs analyses currently.
  if (Parallel::TrajComm().Size() > 1)
    mprintf("Warning: Analysis does not currently use multiple MPI threads.\n");
  if (Parallel::TrajComm().Master())
# endif
    err = analysisList_.DoAnalyses();
  analysis_time_.Stop();
# ifdef MPI
  if (Parallel::World().CheckError( err )) err = 1;
# endif
  mprintf("TIME: Analyses took %.4f seconds.\n", analysis_time_.Total());
  // If all Analyses completed successfully, clean up analyses.
  if ( err == 0) 
    analysisList_.Clear();
  return err;
}

// CpptrajState::AddReference()
int CpptrajState::AddReference( std::string const& fname ) {
  return AddReference( fname, ArgList() );
}

// CpptrajState::AddReference()
/** Add specified file/COORDS set as reference. Reference frames are a unique
  * DataSet - they are set up OUTSIDE data set list.
  */
int CpptrajState::AddReference( std::string const& fname, ArgList const& args ) {
  if (fname.empty()) return 1;
  ArgList argIn = args;
  // 'average' keyword is deprecated
  if ( argIn.hasKey("average") ) {
    mprinterr("Error: 'average' for reference is deprecated. Please use\n"
              "Error:   the 'average' action to create averaged coordinates.\n");
    return 1;
  }
  Topology* refParm = 0;
  DataSet_Coords* CRD = 0;
  if (argIn.hasKey("crdset")) {
    CRD = (DataSet_Coords*)DSL_.FindCoordsSet( fname );
    if (CRD == 0) {
      mprinterr("COORDS set with name %s not found.\n", fname.c_str());
      return 1;
    }
  } else {
    // Get topology file.
    refParm = DSL_.GetTopology( argIn );
    if (refParm == 0) {
      mprinterr("Error: Cannot get topology for reference '%s'\n", fname.c_str());
      return 1;
    }
  }
  std::string tag = argIn.GetStringKey("name");
  // Determine if there is a mask expression for stripping reference. // TODO: Remove?
  std::string maskexpr = argIn.GetMaskNext();
  // Check for tag. FIXME: need to do after SetupTrajRead?
  if (tag.empty()) tag = argIn.getNextTag();
  // Set up reference DataSet from file or COORDS set.
  DataSet_Coords_REF* ref = new DataSet_Coords_REF();
  if (ref==0) return 1;
  if (refParm != 0) {
    if (ref->LoadRefFromFile(fname, tag, *refParm, argIn, refDebug_)) return 1;
  } else { // CRD != 0
    int fnum;
    if (argIn.hasKey("lastframe"))
      fnum = (int)CRD->Size()-1;
    else
      fnum = argIn.getNextInteger(1) - 1;
    mprintf("\tSetting up reference from COORDS set '%s', frame %i\n",
            CRD->legend(), fnum+1);
    if (ref->SetRefFromCoords(CRD, tag, fnum)) return 1;
  }
  // If a mask expression was specified, strip to match the expression.
  if (!maskexpr.empty()) {
    if (ref->StripRef( maskexpr )) return 1;
  }
  // Add DataSet to main DataSetList.
  if (DSL_.AddSet( ref )) return 1; 
  return 0;
}

// CpptrajState::AddTopology()
/** Add specified file(s) as Topology. Topologies are a unique
  * DataSet - they are set up OUTSIDE data set list.
  */
int CpptrajState::AddTopology( std::string const& fnameIn, ArgList const& args ) {
  if (fnameIn.empty()) return 1;
  File::NameArray fnames = File::ExpandToFilenames( fnameIn );
  if (fnames.empty()) {
    mprinterr("Error: '%s' corresponds to no files.\n");
    return 1;
  }
  ArgList argIn = args;
  std::string tag = argIn.GetStringKey("name");
  // Determine if there is a mask expression for stripping. // TODO: Remove?
  std::string maskexpr = argIn.GetMaskNext();
  // Check for tag if 'name' not specified.
  if (tag.empty()) tag = argIn.getNextTag();
  for (File::NameArray::const_iterator fname = fnames.begin(); fname != fnames.end(); ++fname)
  {
    MetaData md(*fname, tag, -1);
    DataSet* set = DSL_.CheckForSet( md );
    if (set != 0) {
      mprintf("Warning: Topology '%s' already present.\n", set->legend());
      mprintf("Warning:   To load the same topology file multiple times use tags,\n"
              "Warning:   e.g. `parm <file> [tag]`.\n");
    } else {
      // Create Topology DataSet
      DataSet_Topology* ds = (DataSet_Topology*)DSL_.AddSet(DataSet::TOPOLOGY, md);
      if (ds == 0) { 
        if (exitOnError_) return 1;
      } else {
        if (ds->LoadTopFromFile(argIn, topDebug_)) {
          DSL_.RemoveSet( ds );
          if (exitOnError_) return 1;
        }
        // If a mask expression was specified, strip to match the expression.
        if (!maskexpr.empty()) {
          if (ds->StripTop( maskexpr )) return 1;
        }
      }
    }
    // TODO: Set active top?
  }
  return 0;
}

int CpptrajState::AddTopology( Topology const& top, std::string const& parmname ) {
  DataSet_Topology* ds = (DataSet_Topology*)DSL_.AddSet(DataSet::TOPOLOGY, parmname);
  if (ds == 0) return 1;
  ds->SetTop( top );
  return 0;
}
