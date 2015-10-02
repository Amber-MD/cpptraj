#include "CpptrajState.h"
#include "CpptrajStdio.h"
#include "MpiRoutines.h" // worldrank
#include "Action_CreateCrd.h" // in case default COORDS need to be created
#include "Timer.h"
#include "DataSet_Coords_REF.h" // AddReference
#include "ProgressBar.h"
#ifdef MPI
# include "DataSet_Coords_TRJ.h"
# ifdef TIMER
#   include "EnsembleIn.h"
# endif
#endif

// CpptrajState::AddTrajin()
int CpptrajState::AddTrajin( ArgList& argIn, bool isEnsemble ) {
  std::string fname = argIn.GetStringNext();
  Topology* tempParm = parmFileList_.GetParm( argIn );
  if (isEnsemble) {
    if ( trajinList_.AddEnsemble( fname, tempParm, argIn ) ) return 1;
  } else {
    if ( trajinList_.AddTrajin( fname, tempParm, argIn ) ) return 1;
  }
  return 0;
}

// CpptrajState::AddTrajin()
int CpptrajState::AddTrajin( std::string const& fname ) {
  if ( trajinList_.AddTrajin( fname, parmFileList_.GetParm(0), ArgList() ) ) return 1;
  return 0;
}

int CpptrajState::AddOutputTrajectory( ArgList& argIn ) {
  std::string fname = argIn.GetStringNext();
  Topology* top = parmFileList_.GetParm( argIn );
  return trajoutList_.AddTrajout( fname, argIn, top );
}

int CpptrajState::AddOutputTrajectory( std::string const& fname ) {
  // Should this use the last Topology instead?
  return trajoutList_.AddTrajout( fname, ArgList(), parmFileList_.GetParm(0) );
}
// -----------------------------------------------------------------------------
int CpptrajState::WorldSize() { return worldsize; }

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
  if ( enabled[L_REF]      ) ReferenceInfo();
  if ( enabled[L_TRAJOUT]  ) trajoutList_.List();
  if ( enabled[L_PARM]     ) parmFileList_.List();
  if ( enabled[L_ANALYSIS] ) analysisList_.List();
  if ( enabled[L_DATAFILE] ) DFL_.List();
  if ( enabled[L_DATASET]  ) DSL_.List();
  return 0;
}

/** Set debug level of specified lists. */
int CpptrajState::SetListDebug( ArgList& argIn ) {
  debug_ = argIn.getNextInteger(0);
  std::vector<bool> enabled = ListsFromArg( argIn, true );
  if ( enabled[L_ACTION]   ) actionList_.SetDebug( debug_ );
  if ( enabled[L_TRAJIN]   ) trajinList_.SetDebug( debug_ );
//  if ( enabled[L_REF]      ) refFrames_.SetDebug( debug_ );
  if ( enabled[L_TRAJOUT]  ) trajoutList_.SetDebug( debug_ );
  if ( enabled[L_PARM]     ) parmFileList_.SetDebug( debug_ );
  if ( enabled[L_ANALYSIS] ) analysisList_.SetDebug( debug_ );
  if ( enabled[L_DATAFILE] ) DFL_.SetDebug( debug_ );
  if ( enabled[L_DATASET]  ) DSL_.SetDebug( debug_ );
  return 0;
}

/** Clear specified lists */
int CpptrajState::ClearList( ArgList& argIn ) {
  std::vector<bool> enabled = ListsFromArg( argIn, false );
  if ( enabled[L_ACTION]   ) actionList_.Clear();
  if ( enabled[L_TRAJIN]   ) trajinList_.Clear();
//  if ( enabled[L_REF]      ) refFrames_.Clear();
  if ( enabled[L_TRAJOUT]  ) trajoutList_.Clear();
  if ( enabled[L_PARM]     ) parmFileList_.Clear();
  if ( enabled[L_ANALYSIS] ) analysisList_.Clear();
  if ( enabled[L_DATAFILE] ) DFL_.Clear();
  if ( enabled[L_DATASET]  ) DSL_.Clear();
  return 0;
}

/** List reference information. */
void CpptrajState::ReferenceInfo() const {
  DSL_.ListReferenceFrames();
  if (activeRef_ != 0)
    mprintf("\tActive reference frame for distance-based masks is '%s'\n", activeRef_->legend());
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
      DFL_.RemoveDataSet( *ds );
      DSL_.RemoveSet( *ds );
    }
  }
  return 0;
}

// CpptrajState::TrajLength()
// NOTE: MMPBSA.py relies on this.
int CpptrajState::TrajLength( std::string const& topname, 
                              std::vector<std::string> const& trajinFiles)
{
  if (parmFileList_.AddParmFile( topname )) return 1;
  for (std::vector<std::string>::const_iterator trajinName = trajinFiles.begin();
                                                trajinName != trajinFiles.end();
                                                ++trajinName)
    if (AddTrajin( *trajinName )) return 1;
  loudPrintf("Frames: %i\n", trajinList_.MaxFrames());
  return 0;
}

// -----------------------------------------------------------------------------
// CpptrajState::Run()
int CpptrajState::Run() {
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
      ArgList crdcmd("createcrd _DEFAULTCRD_");
      crdcmd.MarkArg(0);
      if (AddAction( Action_CreateCrd::Alloc, crdcmd ))
        return 1;
    }
  }
  mprintf("---------- RUN BEGIN -------------------------------------------------\n");
  if (trajinList_.empty()) 
    mprintf("Warning: No input trajectories specified.\n");
  else if (actionList_.Empty() && trajoutList_.Empty())
    mprintf("Warning: No actions/output trajectories specified.\n");
  else {
    switch ( trajinList_.Mode() ) {
#     ifdef MPI
      case TrajinList::NORMAL   : err = RunParallel(); break;
#     else
      case TrajinList::NORMAL   : err = RunNormal(); break;
#     endif
      case TrajinList::ENSEMBLE : err = RunEnsemble(); break;
      case TrajinList::UNDEFINED: break;
    }
    // Clean up Actions if run completed successfully.
    if (err == 0) {
      actionList_.Clear();
      trajoutList_.Clear();
      DSL_.SetDataSetsPending(false);
    }
  }
  // Run Analyses if any are specified.
  if (err == 0)
    err = RunAnalyses();
  DSL_.List();
  // Print DataFile information and write DataFiles
  DFL_.List();
  MasterDataFileWrite();
  mprintf("---------- RUN END ---------------------------------------------------\n");
  return err;
}

// CpptrajState::ActiveReference()
Frame CpptrajState::ActiveReference() const {
  if (activeRef_ == 0)
    return Frame();
  else
    return activeRef_->RefFrame();
}

// -----------------------------------------------------------------------------
// CpptrajState::RunEnsemble()
int CpptrajState::RunEnsemble() {
  Timer init_time;
  init_time.Start();
  FrameArray FrameEnsemble;
  FramePtrArray SortedFrames;

  mprintf("\nINPUT ENSEMBLE:\n");
  // Ensure all ensembles are of the same size
  int ensembleSize = -1;
  for (TrajinList::ensemble_it traj = trajinList_.ensemble_begin(); 
                               traj != trajinList_.ensemble_end(); ++traj) 
  {
    if (ensembleSize == -1) {
      ensembleSize = (*traj)->EnsembleCoordInfo().EnsembleSize();
#     ifdef MPI
      // TODO: Eventually try to divide ensemble among MPI threads?
      if (worldsize != ensembleSize) {
        mprinterr("Error: Ensemble size (%i) does not match # of MPI threads (%i).\n",
                  ensembleSize, worldsize);
        return 1;
      }
#     endif
    } else if (ensembleSize != (*traj)->EnsembleCoordInfo().EnsembleSize()) {
      mprinterr("Error: Ensemble size (%i) does not match first ensemble size (%i).\n",
                (*traj)->EnsembleCoordInfo().EnsembleSize(), ensembleSize);
      return 1;
    }
  }
  mprintf("  Ensemble size is %i\n", ensembleSize);
  // Allocate space to hold position of each incoming frame in replica space.
# ifdef MPI
  // Only two frames needed; one for reading, one for receiving.
  SortedFrames.resize( 2 );
  FrameEnsemble.resize( 2 );
# else
  SortedFrames.resize( ensembleSize );
  FrameEnsemble.resize( ensembleSize );
# endif
  //FrameEnsemble.SetupFrames( TrajParm()->Atoms(), cInfo_ );

  // At this point all ensembles should match (i.e. same map etc.)
  trajinList_.FirstEnsembleReplicaInfo();

  // Calculate frame division among trajectories
  trajinList_.List();
  // Parameter file information
  parmFileList_.List();
  // Print reference information 
  ReferenceInfo();
  // Use separate TrajoutList. Existing trajout in current TrajoutList
  // will be converted to ensemble trajout.
  EnsembleOutList ensembleOut;
  if (!trajoutList_.Empty()) {
    if (trajoutList_.MakeEnsembleTrajout(ensembleOut, ensembleSize))
      return 1;
    mprintf("\nENSEMBLE OUTPUT TRAJECTORIES (Numerical filename"
            " suffix corresponds to above map):\n");
    parallel_barrier();
    ensembleOut.List();
    parallel_barrier();
  }
  // Allocate DataSets in the master DataSetList based on # frames to be read
  DSL_.AllocateSets( trajinList_.MaxFrames() );
# ifdef MPI
  // Each thread will process one member of the ensemble, so local ensemble
  // size is effectively 1.
  ensembleSize = 1;
# endif
  // Allocate an ActionList for each member of the ensemble.
  std::vector<ActionList*> ActionEnsemble( ensembleSize );
  ActionEnsemble[0] = &actionList_;
  for (int member = 1; member < ensembleSize; member++)
    ActionEnsemble[member] = new ActionList();
  // If we are on a single thread, give each member its own copy of the
  // current topology address. This way if topology is modified by a member,
  // e.g. in strip or closest, subsequent members wont be trying to modify 
  // an already-modified topology.
  std::vector<Topology*> EnsembleParm( ensembleSize );
  // Give each member its own copy of current frame address. This way if 
  // frame is modified by a member things like trajout know about it.
  FramePtrArray CurrentFrames( ensembleSize );
# ifdef MPI
  // Make all sets not in an ensemble a member of this thread.
  DSL_.MakeDataSetsEnsemble( worldrank );
  // This tells all DataFiles to append member number.
  DFL_.MakeDataFilesEnsemble( worldrank );
  // Actions have already been set up for this ensemble.
# else
  // Make all sets not in an ensemble part of member 0.
  DSL_.MakeDataSetsEnsemble( 0 );
  // Silence action output for members > 0.
  if (debug_ == 0) SetWorldSilent( true ); 
  // Set up Actions for each ensemble member > 0.
  for (int member = 1; member < ensembleSize; ++member) {
    // All DataSets that will be set up will be part of this ensemble 
    DSL_.SetEnsembleNum( member );
    // Initialize actions for this ensemble member based on original actionList_
    if (!actionList_.Empty()) {
      if (debug_ > 0) mprintf("***** ACTIONS FOR ENSEMBLE MEMBER %i:\n", member);
      for (int iaction = 0; iaction < actionList_.Naction(); iaction++) { 
        // Create new arg list from original command string.
        ArgList command( actionList_.CmdString(iaction) );
        command.MarkArg(0); // TODO: Create separate CommandArg class?
        // Attempt to add same action to this ensemble. 
        if (ActionEnsemble[member]->AddAction( actionList_.ActionAlloc(iaction), 
                                                command, &parmFileList_, &DSL_, &DFL_ ))
            return 1;
      }
    }
  }
  SetWorldSilent( false );
# endif
  init_time.Stop();
  // Re-enable output
  SetWorldSilent( false );
  mprintf("TIME: Run Initialization took %.4f seconds.\n", init_time.Total()); 
  // ========== A C T I O N  P H A S E ==========
  int lastPindex=-1;          // Index of the last loaded parm file
  int readSets = 0;
  int actionSet = 0;
  bool hasVelocity = false;
# ifdef TIMER
  Timer trajin_time;
  Timer setup_time;
  Timer actions_time;
  Timer trajout_time;
# endif
  Timer frames_time;
  frames_time.Start();
  // Loop over every trajectory in trajFileList
  mprintf("\nBEGIN ENSEMBLE PROCESSING:\n");
  ProgressBar progress;
  if (showProgress_)
    progress.SetupProgress( trajinList_.MaxFrames() );
  for ( TrajinList::ensemble_it traj = trajinList_.ensemble_begin();
                                traj != trajinList_.ensemble_end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail
    if ( (*traj)->BeginEnsemble() ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->Traj().Filename().full());
      break;
    }
    // Set current parm from current traj.
    Topology* CurrentParm = (*traj)->Traj().Parm();
    for (int member = 0; member < ensembleSize; ++member)
      EnsembleParm[member] = CurrentParm;
    CoordinateInfo const& currentCoordInfo = (*traj)->EnsembleCoordInfo();
    CurrentParm->SetParmCoordInfo( currentCoordInfo );
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != CurrentParm->Pindex());
#   ifdef TIMER
    setup_time.Start();
#   endif
    // If Parm has changed or trajectory velocity status has changed,
    // reset the frame.
    if (parmHasChanged || (hasVelocity != currentCoordInfo.HasVel()))
      FrameEnsemble.SetupFrames(CurrentParm->Atoms(), currentCoordInfo);
    hasVelocity = currentCoordInfo.HasVel();

    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set active reference for this parm
      CurrentParm->SetReferenceCoords( ActiveReference() );
      // Set up actions for this parm
      bool setupOK = true;
      for (int member = 0; member < ensembleSize; ++member) {
        // Silence action output for all beyond first member.
        if (member > 0)
          SetWorldSilent( true );
        if (ActionEnsemble[member]->SetupActions( &(EnsembleParm[member]) )) {
#         ifdef MPI
          rprintf("Warning: Ensemble member %i: Could not set up actions for %s: skipping.\n",
                  worldrank,EnsembleParm[member]->c_str());
#         else
          mprintf("Warning: Ensemble member %i: Could not set up actions for %s: skipping.\n",
                  member,EnsembleParm[member]->c_str());
#         endif
          setupOK = false;
        }
      }
      // Re-enable output
      SetWorldSilent( false );
      if (!setupOK) continue;
      // Set up any related output trajectories.
      // TODO: Currently assuming topology is always modified the same
      //       way for all actions. If this behavior ever changes the
      //       following line will cause undesireable behavior.
      ensembleOut.SetupEnsembleOut( EnsembleParm[0] );
      lastPindex = CurrentParm->Pindex();
    }
#   ifdef TIMER
    setup_time.Stop();
#   endif
    // Loop over every collection of frames in the ensemble
    (*traj)->Traj().PrintInfoLine();
#   ifdef TIMER
    trajin_time.Start();
    bool readMoreFrames = (*traj)->GetNextEnsemble(FrameEnsemble, SortedFrames);
    trajin_time.Stop();
    while ( readMoreFrames )
#   else
    while ( (*traj)->GetNextEnsemble(FrameEnsemble, SortedFrames) )
#   endif
    {
      if (!(*traj)->BadEnsemble()) {
        bool suppress_output = false;
        for (int member = 0; member != ensembleSize; ++member) {
          // Since Frame can be modified by actions, save original and use CurrentFrame
          Frame* CurrentFrame = SortedFrames[member];
          //rprintf("DEBUG: CurrentFrame=%x SortedFrames[0]=%x\n",CurrentFrame, SortedFrames[0]);
          if ( CurrentFrame->CheckCoordsInvalid() )
            rprintf("Warning: Ensemble member %i frame %i may be corrupt.\n",
                    member, (*traj)->Traj().Counter().PreviousFrameNumber()+1);
#         ifdef TIMER
          actions_time.Start();
#         endif
          // Perform Actions on Frame
          suppress_output = ActionEnsemble[member]->DoActions(&CurrentFrame, actionSet);
          CurrentFrames[member] = CurrentFrame;
#         ifdef TIMER
          actions_time.Stop();
#         endif
        } // END loop over actions
        // Do Output
        if (!suppress_output) {
#         ifdef TIMER
          trajout_time.Start();
#         endif 
          if (ensembleOut.WriteEnsembleOut(actionSet, CurrentFrames))
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
      readMoreFrames = (*traj)->GetNextEnsemble(FrameEnsemble, SortedFrames);
      trajin_time.Stop();
#     endif
    }

    // Close the trajectory file
    (*traj)->EndEnsemble();
    // Update how many frames have been processed.
    readSets += (*traj)->Traj().Counter().NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  mprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);
  frames_time.Stop();
  frames_time.WriteTiming(0," Trajectory processing:");
  mprintf("TIME: Avg. throughput= %.4f frames / second.\n",
          (double)readSets / frames_time.Total());
# ifdef TIMER
  trajin_time.WriteTiming(1,  "Trajectory read:        ", frames_time.Total());
  setup_time.WriteTiming(1,   "Action setup:           ", frames_time.Total());
  actions_time.WriteTiming(1, "Action frame processing:", frames_time.Total());
  trajout_time.WriteTiming(1, "Trajectory output:      ", frames_time.Total());
# ifdef MPI
  Ensemble::TimingData(trajin_time.Total());
# endif
# endif

  // Close output trajectories
  ensembleOut.CloseEnsembleOut();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  mprintf("\nENSEMBLE ACTION OUTPUT:\n");
  for (int member = 0; member < ensembleSize; ++member)
    ActionEnsemble[member]->Print( );
# ifdef MPI
  // Sync DataSets across all threads. 
  //DSL_.SynchronizeData(); // NOTE: Disabled, trajs are not currently divided.
# endif
  // Clean up ensemble action lists
  for (int member = 1; member < ensembleSize; member++)
    delete ActionEnsemble[member];

  return 0;
}
#ifdef MPI
// -----------------------------------------------------------------------------
/** Process trajectories in trajinList in parallel. */
int CpptrajState::RunParallel() {
  // Print information.
  parmFileList_.List();
  trajinList_.List();
  ReferenceInfo();

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
  if (parallel_check_error( err )) {
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
  if (parallel_check_error( err )) return 1;

  // Divide frames among threads.
  int my_start = 0;
  int my_stop = 0;
//  int my_frames = 0;
  int maxFrames = (int)input_traj.Size();
  int frames_per_thread = maxFrames / worldsize;
  int remainder         = maxFrames % worldsize;
  int my_frames         = frames_per_thread + (int)(worldrank < remainder);
  // Figure out where this thread starts and stops
  for (int rank = 0; rank != worldrank; rank++)
    if (rank < remainder)
      my_start += (frames_per_thread + 1);
    else
      my_start += (frames_per_thread);
  my_stop = my_start + my_frames;
  // Store how many frames each rank will process (only needed by master).
  std::vector<int> rank_frames( worldsize, frames_per_thread );
  for (int i = 0; i != worldsize; i++)
    if (i < remainder) rank_frames[i]++;
  if (worldrank == 0) {
    for (int i = 0; i != worldsize; i++)
      mprintf("Thread %i will process %i frames.\n", i, rank_frames[i]);
  }
/*  if (worldrank == 0) {
    mprintf("Max frames is %zu\n", input_traj.Size());
    int maxFrames = (int)input_traj.Size();
    int frames_per_thread = maxFrames / worldsize;
    int remainder         = maxFrames % worldsize;
    int* rank_frames = new int[ worldsize * 3 ];
    std::fill( rank_frames, rank_frames + worldsize, frames_per_thread );
    int* rank_start = rank_frames + worldsize;
    int* rank_stop  = rank_start  + worldsize;
    std::fill( rank_start, rank_start + (2*worldsize), 0 ); // TODO Necessary?
    for (int rank = 0; rank != worldsize; rank++) {
      if (rank > 0)
        rank_start[rank] = rank_stop[rank-1];
      if (rank < remainder)
        rank_frames[rank]++;
      rank_stop[rank] = rank_start[rank] + rank_frames[rank];
      //rprintf("Start %i Stop %i\n", rank_start[rank], rank_stop[rank]);
      if (rank > 0) {
        parallel_send( rank_start  + rank, 1, PARA_INT, rank, 1000 );
        parallel_send( rank_stop   + rank, 1, PARA_INT, rank, 1001 );
        parallel_send( rank_frames + rank, 1, PARA_INT, rank, 1002 );
      }
    }
    my_stop = rank_stop[0];
    my_frames = rank_frames[0];
    delete[] rank_frames;
  } else {
    // Get start, stop, and total frames from master
    parallel_recv( &my_start,  1, PARA_INT, 0, 1000 );
    parallel_recv( &my_stop,   1, PARA_INT, 0, 1001 );
    parallel_recv( &my_frames, 1, PARA_INT, 0, 1002 );
  }*/
  rprintf("Start %i Stop %i Frames %i\n", my_start, my_stop, my_frames);
  parallel_barrier();

  // Disable trajectory output for now. FIXME
  if (!trajoutList_.Empty()) {
    mprinterr("Error: trajout not yet supported in parallel.\n");
    return 1;
  }

  // Allocate DataSets in DataSetList based on # frames read by this thread.
  DSL_.AllocateSets( my_frames );

  // ----- SETUP PHASE ---------------------------
  CoordinateInfo const& currentCoordInfo = input_traj.CoordsInfo();
  Topology* CurrentParm = (Topology*)&(input_traj.Top()); // TODO fix cast
  CurrentParm->SetParmCoordInfo( currentCoordInfo );
  CurrentParm->SetReferenceCoords( ActiveReference() );
  err = actionList_.SetupActions( &CurrentParm );
  if (parallel_check_error( err )) {
    mprinterr("Error: Could not set up actions for '%s'\n", CurrentParm->c_str());
    return 1;
  }
  // TODO Set up output trajectories

  // ----- ACTION PHASE --------------------------
  ProgressBar progress;
  if (showProgress_)
    progress.SetupProgress( my_frames );
  Frame TrajFrame = input_traj.AllocateFrame();
  int actionSet = 0; // Internal data frame
  for (int set = my_start; set != my_stop; set++, actionSet++) {
    input_traj.GetFrame(set, TrajFrame);
    if (TrajFrame.CheckCoordsInvalid()) // TODO actual frame #
      rprintf("Warning: Set %i coords 1 & 2 overlap at origin; may be corrupt.\n", set + 1);
    Frame* CurrentFrame = &TrajFrame;
    //bool suppress_output =
    actionList_.DoActions(&CurrentFrame, actionSet);
    // TODO trajout stuff
    if (showProgress_) progress.Update( actionSet );
  }
  parallel_barrier();
  // Sync data sets to master thread
  Timer time_sync;
  time_sync.Start();
  DSL_.SynchronizeData( input_traj.Size(), rank_frames );
  time_sync.Stop();
  time_sync.WriteTiming(1, "Data set sync");
  mprintf("\nACTION OUTPUT:\n");
  // Only call print for master
  if (worldrank == 0)
    actionList_.Print( );
  parallel_barrier();
  return 0;
}
#endif
// -----------------------------------------------------------------------------
// CpptrajState::RunNormal()
/** Process trajectories in trajinList. Each frame in trajinList is sent
 *  to the actions in actionList for processing.
 */
int CpptrajState::RunNormal() {
  int actionSet=0;            // Internal data frame
  int readSets=0;             // Number of frames actually read
  int lastPindex=-1;          // Index of the last loaded parm file
  Frame TrajFrame;            // Original Frame read in from traj

  // ========== S E T U P   P H A S E ========== 
  Timer init_time;
  init_time.Start();
  // Parameter file information
  parmFileList_.List();
  // Input coordinate file information
  trajinList_.List();
  // Print reference information
  ReferenceInfo(); 
  // Output traj
  trajoutList_.List();
  // Allocate DataSets in the master DataSetList based on # frames to be read
  DSL_.AllocateSets( trajinList_.MaxFrames() );
  init_time.Stop();
  mprintf("TIME: Run Initialization took %.4f seconds.\n", init_time.Total());
  
  // ========== A C T I O N  P H A S E ==========
  // Loop over every trajectory in trajFileList
# ifdef TIMER
  Timer trajin_time;
  Timer setup_time;
  Timer actions_time;
  Timer trajout_time;
# endif
  Timer frames_time;
  frames_time.Start();
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
    CoordinateInfo const& currentCoordInfo = (*traj)->TrajCoordInfo();
    Topology* CurrentParm = (*traj)->Traj().Parm();
    CurrentParm->SetParmCoordInfo( currentCoordInfo );
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != CurrentParm->Pindex());
#   ifdef TIMER
    setup_time.Start();
#   endif
    // If Parm has changed or trajectory frame has changed, reset the frame.
    if (parmHasChanged || 
        (TrajFrame.HasVelocity() != currentCoordInfo.HasVel()) ||
        ((int)TrajFrame.RemdIndices().size() != currentCoordInfo.ReplicaDimensions().Ndims()))
      TrajFrame.SetupFrameV(CurrentParm->Atoms(), currentCoordInfo);
    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set active reference for this parm
      CurrentParm->SetReferenceCoords( ActiveReference() );
      // Set up actions for this parm
      if (actionList_.SetupActions( &CurrentParm )) {
        mprintf("WARNING: Could not set up actions for %s: skipping.\n",
                CurrentParm->c_str());
        continue;
      }
      // Set up any related output trajectories 
      trajoutList_.SetupTrajout( CurrentParm );
      lastPindex = CurrentParm->Pindex();
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
      // Check that coords are valid.
      if ( TrajFrame.CheckCoordsInvalid() )
        mprintf("Warning: Frame %i coords 1 & 2 overlap at origin; may be corrupt.\n",
                (*traj)->Traj().Counter().PreviousFrameNumber()+1);
        // Since Frame can be modified by actions, save original and use CurrentFrame
        Frame* CurrentFrame = &TrajFrame;
        // Perform Actions on Frame
#       ifdef TIMER
        actions_time.Start();
#       endif
        bool suppress_output = actionList_.DoActions(&CurrentFrame, actionSet);
#       ifdef TIMER
        actions_time.Stop();
#       endif
        // Do Output
        if (!suppress_output) {
#         ifdef TIMER
          trajout_time.Start();
#         endif
          if (trajoutList_.WriteTrajout(actionSet, *CurrentFrame)) {
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
  frames_time.Stop();
  frames_time.WriteTiming(0," Trajectory processing:");
  mprintf("TIME: Avg. throughput= %.4f frames / second.\n", 
          (double)readSets / frames_time.Total());
# ifdef TIMER
  trajin_time.WriteTiming(1,  "Trajectory read:        ", frames_time.Total());
  setup_time.WriteTiming(1,   "Action setup:           ", frames_time.Total());
  actions_time.WriteTiming(1, "Action frame processing:", frames_time.Total());
  trajout_time.WriteTiming(1, "Trajectory output:      ", frames_time.Total());
# endif
  // Close output trajectories.
  trajoutList_.CloseTrajout();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  mprintf("\nACTION OUTPUT:\n");
  actionList_.Print( );
# ifdef MPI
  // Sync DataSets across all threads. 
  //DSL_.SynchronizeData(); // NOTE: Disabled, trajs are not currently divided.
# endif

  return 0;
}

// -----------------------------------------------------------------------------
// CpptrajState::MasterDataFileWrite()
// FIXME: If MPI ever used for trajin mode this may have to be protected.
/** Trigger write of all pending DataFiles. When in parallel ensemble mode,
  * each member of the ensemble will write data to separate files with 
  * numeric extensions.
  */
void CpptrajState::MasterDataFileWrite() { DFL_.WriteAllDF(); }

// CpptrajState::RunAnalyses()
int CpptrajState::RunAnalyses() {
  if (analysisList_.Empty()) return 0;
  Timer analysis_time;
  analysis_time.Start();
  // Only master performs analyses currently.
  int err = 0;
  if (worldrank == 0)
    err = analysisList_.DoAnalyses();
  analysis_time.Stop();
# ifdef MPI
  if (parallel_check_error( err )) err = 1;
# endif
  mprintf("TIME: Analyses took %.4f seconds.\n", analysis_time.Total());
  // If all Analyses completed successfully, clean up analyses.
  if ( err == 0) 
    analysisList_.Clear();
  return err;
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
    refParm = parmFileList_.GetParm( argIn );
    if (refParm == 0) {
      mprinterr("Error: Cannot get topology for reference '%s'\n", fname.c_str());
      return 1;
    }
  }
  // Determine if there is a mask expression for stripping reference. // TODO: Remove?
  std::string maskexpr = argIn.GetMaskNext();
  // Check for tag. FIXME: need to do after SetupTrajRead?
  std::string tag = argIn.getNextTag();
  // Set up reference DataSet from file or COORDS set.
  DataSet_Coords_REF* ref = new DataSet_Coords_REF();
  if (ref==0) return 1;
  if (refParm != 0) {
    if (ref->LoadRefFromFile(fname, tag, *refParm, argIn, debug_)) return 1;
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
  // Set default reference if not already set.
  if (activeRef_ == 0) activeRef_ = ref;
  return 0;
}
