#include "CpptrajState.h"
#include "CpptrajStdio.h"
#include "FrameArray.h" // for ensemble
#include "Trajin_Multi.h" // for ensemble
#include "MpiRoutines.h" // worldrank
#include "Action_CreateCrd.h" // in case default COORDS need to be created
#ifdef TIMER
# include "Timer.h"
#endif

int CpptrajState::WorldSize() { return worldsize; }

/** Select lists from ArgList */
std::vector<bool> CpptrajState::ListsFromArg( ArgList& argIn, bool allowEnableAll ) {
  std::vector<bool> enabled( (int)N_LISTS );
  enabled[L_ACTION]   = argIn.hasKey("actions");
  enabled[L_TRAJIN]   = argIn.hasKey("trajin");
  enabled[L_REF]      = argIn.hasKey("ref");
  enabled[L_TRAJOUT]  = argIn.hasKey("trajout");
  enabled[L_PARM]     = argIn.hasKey("parm");
  enabled[L_ANALYSIS] = argIn.hasKey("analysis");
  enabled[L_DATAFILE] = argIn.hasKey("datafile");
  enabled[L_DATASET]  = argIn.hasKey("dataset");
  if (!allowEnableAll) return enabled;
  // If nothing is enabled, set all enabled
  bool nothing_enabled = true;
  for (std::vector<bool>::iterator en = enabled.begin(); en != enabled.end(); ++en)
    if (*en) {
      nothing_enabled = false;
      break;
    }
  if (nothing_enabled) enabled.assign( (int)N_LISTS, true );
  return enabled;
}

/** List all members of specified lists */
int CpptrajState::ListAll( ArgList& argIn ) {
  std::vector<bool> enabled = ListsFromArg( argIn, true );
  if ( enabled[L_ACTION]   ) actionList_.List();
  if ( enabled[L_TRAJIN]   ) trajinList_.List();
  if ( enabled[L_REF]      ) refFrames_.List();//{mprintf("\nREFERENCE COORDS:\n");refFrames_.List();}
  if ( enabled[L_TRAJOUT]  ) {mprintf("\nOUTPUT TRAJECTORIES:\n");trajoutList_.List();}
  if ( enabled[L_PARM]     ) parmFileList_.List();
  if ( enabled[L_ANALYSIS] ) analysisList_.List();
  if ( enabled[L_DATAFILE] ) DFL_.List();
  if ( enabled[L_DATASET]  ) {mprintf("\nDATASETS:\n");DSL_.List();}
  return 0;
}

/** Set debug level of specified lists. */
int CpptrajState::SetListDebug( ArgList& argIn ) {
  std::vector<bool> enabled = ListsFromArg( argIn, true );
  debug_ = argIn.getNextInteger(0);
  if ( enabled[L_ACTION]   ) actionList_.SetDebug( debug_ );
  if ( enabled[L_TRAJIN]   ) trajinList_.SetDebug( debug_ );
  if ( enabled[L_REF]      ) refFrames_.SetDebug( debug_ );
  if ( enabled[L_TRAJOUT]  ) trajoutList_.SetDebug( debug_ );
  if ( enabled[L_PARM]     ) parmFileList_.SetDebug( debug_ );
  if ( enabled[L_ANALYSIS] ) analysisList_.SetDebug( debug_ );
  if ( enabled[L_DATAFILE] ) DFL_.SetDebug( debug_ );
  if ( enabled[L_DATASET]  ) DSL_.SetDebug( debug_ );
  return 0;
}

/** Clear specified lists */
int CpptrajState::ClearList( ArgList& argIn ) {
  std::vector<bool> enabled = ListsFromArg( argIn, argIn.hasKey("all") );
  if ( enabled[L_ACTION]   ) actionList_.Clear();
  if ( enabled[L_TRAJIN]   ) trajinList_.Clear();
  if ( enabled[L_REF]      ) refFrames_.Clear();
  if ( enabled[L_TRAJOUT]  ) trajoutList_.Clear();
  if ( enabled[L_PARM]     ) parmFileList_.Clear();
  if ( enabled[L_ANALYSIS] ) analysisList_.Clear();
  if ( enabled[L_DATAFILE] ) DFL_.Clear();
  if ( enabled[L_DATASET]  ) DSL_.Clear();
  return 0;
}

// CpptrajState::MaskString()
int CpptrajState::MaskString( std::string const& maskexpr ) {
  Topology* parm = parmFileList_.GetParm( 0 );
  if (parm == 0) {
    mprinterr("Error: No topologies loaded.\n");
    return 1;
  }
  parm->PrintAtomInfo( maskexpr );
  return 0;
}

// -----------------------------------------------------------------------------
// CpptrajState::Run()
int CpptrajState::Run() {
# ifdef TIMER
  Timer total_time;
  total_time.Start();
# endif
  ++nrun_;
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
      ArgList crdcmd("createcrd _DEFAULTCRD_");
      crdcmd.MarkArg(0);
      if (AddAction( Action_CreateCrd::Alloc, crdcmd ))
        return 1;
    }
  }

  int err = 0;
  switch ( trajinList_.Mode() ) {
    case TrajinList::NORMAL   :
      err = RunNormal();
      break;
    case TrajinList::ENSEMBLE :
      // No Analysis will be run. Warn user if analyses are defined.
      if (!analysisList_.Empty())
        mprintf("Warning: In ensemble mode, Analysis will not be performed.\n");
      err = RunEnsemble(); 
      break;
    default:
      // No trajectories loaded; If analyses are defined, try to run them.
      if (!analysisList_.Empty()) {
        RunAnalyses();
        MasterDataFileWrite();
      }
  }
  // Clean up Actions.
  actionList_.Clear();
# ifdef TIMER
  total_time.Stop();
  mprintf("TIME: Total Run execution time: %.4f seconds.\n", total_time.Total());
# endif
  return err;
}

// CpptrajState::RunEnsemble()
int CpptrajState::RunEnsemble() {
  FrameArray FrameEnsemble;

  mprintf("\nINPUT ENSEMBLE:\n");
  // Ensure all ensembles are of the same size
  int ensembleSize = -1;
  for (TrajinList::const_iterator traj = trajinList_.begin(); traj != trajinList_.end(); ++traj) 
  {
    Trajin_Multi* mtraj = (Trajin_Multi*)*traj;
    if (ensembleSize == -1) {
      ensembleSize = mtraj->EnsembleSize();
#     ifdef MPI
      // TODO: Eventually try to divide ensemble among MPI threads?
      if (worldsize != ensembleSize) {
        mprinterr("Error: Ensemble size (%i) does not match # of MPI threads (%i).\n",
                  ensembleSize, worldsize);
        return 1;
      }
#     endif
    } else if (ensembleSize != mtraj->EnsembleSize()) {
      mprinterr("Error: Ensemble size (%i) does not match first ensemble size (%i).\n",
                mtraj->EnsembleSize(), ensembleSize);
      return 1;
    }
    // Perform ensemble setup - this also resizes FrameEnsemble
    if ( mtraj->EnsembleSetup( FrameEnsemble ) ) return 1;
  }
  mprintf("  Ensemble size is %i\n", ensembleSize); 
  // At this point all ensembles should match (i.e. same map etc.)
  ((Trajin_Multi*)(trajinList_.front()))->EnsembleInfo();

  // Calculate frame division among trajectories
  trajinList_.List();
  int maxFrames = trajinList_.MaxFrames();
  // Parameter file information
  parmFileList_.List();
  // Print reference information 
  refFrames_.List();
# ifdef MPI
  // Each thread will process one member of the ensemble, so total ensemble
  // size is effectively 1.
  ensembleSize = 1;
# endif
  // Allocate an ActionList, TrajoutList, and DataSetList for each
  // member of the ensemble. Use separate DataFileList.
  std::vector<ActionList> ActionEnsemble( ensembleSize );
  std::vector<TrajoutList> TrajoutEnsemble( ensembleSize );
  std::vector<DataSetList> DataSetEnsemble( ensembleSize );
  DataFileList DataFileEnsemble;
# ifdef MPI
  DataFileEnsemble.SetEnsembleMode( worldrank );
# endif

  // Set up output trajectories for each member of the ensemble
  for (ArgsArray::iterator targ = trajoutArgs_.begin(); targ != trajoutArgs_.end(); ++targ)
  {
#   ifdef MPI
    TrajoutEnsemble[0].AddEnsembleTrajout( *targ, parmFileList_, worldrank );
#   else
    for (int member = 0; member < ensembleSize; ++member) 
      TrajoutEnsemble[member].AddEnsembleTrajout( *targ, parmFileList_, member );
#   endif
  }
  mprintf("\nENSEMBLE OUTPUT TRAJECTORIES (Numerical filename suffix corresponds to above map):\n");
  TrajoutEnsemble[0].List();
  if (debug_ > 0) {
    for (int member = 1; member < ensembleSize; ++member) {
      mprintf("OUTPUT TRAJECTORIES Member %i:\n", member);
      TrajoutEnsemble[member].List();
    }
  }

  // TODO: One loop over member?
  for (int member = 0; member < ensembleSize; ++member) {
    // Set max frames in the data set list and allocate
    DataSetEnsemble[member].SetMax( maxFrames );
    DataSetEnsemble[member].AllocateSets();
    // Initialize actions for this ensemble member based on original actionList_
    if (!actionList_.Empty()) {
      mprintf("***** ACTIONS FOR ENSEMBLE MEMBER %i:\n", member);
      for (int iaction = 0; iaction < actionList_.Naction(); iaction++) { 
        // Create new arg list from original command string.
        ArgList command( actionList_.CmdString(iaction) );
        command.MarkArg(0); // TODO: Create separate CommandArg class?
        // Attempt to add same action to this ensemble. 
        if (ActionEnsemble[member].AddAction( actionList_.ActionAlloc(iaction), 
                                              command, &parmFileList_, &refFrames_,
                                              &(DataSetEnsemble[member]), 
                                              &DataFileEnsemble ))
            return 1;
      }
    }
  }
      
  // ========== A C T I O N  P H A S E ==========
  int lastPindex=-1;          // Index of the last loaded parm file
  int pos = 0;                // Where member should be processed by actions
  int readSets = 0;
  int actionSet = 0;
  bool hasVelocity = false;
  // Loop over every trajectory in trajFileList
  mprintf("\nBEGIN ENSEMBLE PROCESSING:\n");
  for ( TrajinList::const_iterator traj = trajinList_.begin();
                                   traj != trajinList_.end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail 
    if ( (*traj)->BeginTraj(showProgress_) ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->TrajFilename().full());
      break;
    }
    // Set current parm from current traj.
    Topology* CurrentParm = (*traj)->TrajParm();
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != CurrentParm->Pindex());

    // If Parm has changed or trajectory velocity status has changed,
    // reset the frame.
    if (parmHasChanged || (hasVelocity != (*traj)->HasVelocity()))
      FrameEnsemble.SetupFrames(CurrentParm->Atoms(), (*traj)->HasVelocity(),
                                (*traj)->NreplicaDimension());
    hasVelocity = (*traj)->HasVelocity();

    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set active reference for this parm
      CurrentParm->SetReferenceCoords( refFrames_.ActiveReference() );
      // Set up actions for this parm
      bool setupOK = true;
      for (int member = 0; member < ensembleSize; ++member) {
        if (ActionEnsemble[member].SetupActions( &CurrentParm )) {
#         ifdef MPI
          rprintf("Warning: Ensemble member %i: Could not set up actions for %s: skipping.\n",
                  worldrank,CurrentParm->c_str());
#         else
          mprintf("Warning: Ensemble member %i: Could not set up actions for %s: skipping.\n",
                  member,CurrentParm->c_str());
#         endif
          setupOK = false;
        }
      }
      if (!setupOK) continue;
      lastPindex = CurrentParm->Pindex();
    }

    // Loop over every collection of frames in the ensemble
    (*traj)->PrintInfoLine();
    Trajin_Multi* mtraj = (Trajin_Multi*)*traj;
    while ( mtraj->GetNextEnsemble(FrameEnsemble) ) {
      if (!mtraj->BadEnsemble()) {
#       ifdef MPI
        // For MPI, each thread has one ensemble frame. member is 1 if coords
        // had to be sorted, 0 otherwise. pos is always 0.
        int member = mtraj->EnsembleFrameNum();
        pos = 0;
#       else
        // Loop over all members of the ensemble
        for (int member = 0; member < ensembleSize; ++member) {
          // Get this members current position
          pos = mtraj->EnsemblePosition( member );
#       endif
          // Since Frame can be modified by actions, save original and use CurrentFrame
          Frame* CurrentFrame = &(FrameEnsemble[member]);
          // Perform Actions on Frame
          bool suppress_output = ActionEnsemble[pos].DoActions(&CurrentFrame, actionSet);
          // Do Output
          if (!suppress_output) 
            TrajoutEnsemble[pos].Write(actionSet, CurrentParm, CurrentFrame);
#       ifndef MPI
        } // END loop over ensemble
#       endif
      } else {
#       ifdef MPI
        rprinterr("Error: Could not read frame %i for ensemble.\n", actionSet + 1);
#       else
        mprinterr("Error: Could not read frame %i for ensemble.\n", actionSet + 1);
#       endif
      }
      // Increment frame counter
      ++actionSet;
    }

    // Close the trajectory file
    (*traj)->EndTraj();
    // Update how many frames have been processed.
    readSets += (*traj)->NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  rprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);

  // Close output trajectories
  for (int member = 0; member < ensembleSize; ++member)
    TrajoutEnsemble[member].Close();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  mprintf("\nENSEMBLE ACTION OUTPUT:\n");
  for (int member = 0; member < ensembleSize; ++member)
    ActionEnsemble[member].Print( );

  // Sort DataSets and print DataSet information
  // TODO - Also have datafilelist call a sync??
  unsigned int total_data_sets = DataSetEnsemble[0].size();
  mprintf("\nENSEMBLE DATASETS: Each member has %u sets total.\n", total_data_sets);
  for (int member = 0; member < ensembleSize; ++member) {
    //DataSetEnsemble[member].Sync(); // SYNC only necessary when splitting up data
    DataSetEnsemble[member].sort();
    if (total_data_sets != DataSetEnsemble[member].size())
      mprintf("Warning: Ensemble member %i # data sets (%i) does not match member 0 (%i)\n",
              member, DataSetEnsemble[member].size(), total_data_sets);
    if (debug_ > 0)
      DataSetEnsemble[member].List();
  }

  // Print Datafile information
  DataFileEnsemble.List();
  // Print DataFiles. When in parallel ensemble mode, each member of the 
  // ensemble will write data to separate files with numeric extensions. 
  DataFileEnsemble.WriteAllDF();

  return 0;
}

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
# ifdef TIMER
  Timer init_time;
  init_time.Start();
# endif
  // Parameter file information
  parmFileList_.List();
  // Input coordinate file information
  trajinList_.List();
  // Print reference information 
  refFrames_.List();
  // Output traj
  mprintf("\nOUTPUT TRAJECTORIES:\n");
  trajoutList_.List();
  // Allocate DataSets in the master DataSetList based on # frames to be read
  DSL_.AllocateSets();
# ifdef TIMER
  init_time.Stop();
  mprintf("TIME: Run Initialization took %.4f seconds.\n", init_time.Total());
# endif
  
  // ========== A C T I O N  P H A S E ==========
  // Loop over every trajectory in trajFileList
# ifdef TIMER
  Timer actions_time;
  Timer setup_time;
  Timer trajin_time;
  Timer trajout_time;
  Timer frames_time;
  frames_time.Start();
# endif
  mprintf("\nBEGIN TRAJECTORY PROCESSING:\n");
  for ( TrajinList::const_iterator traj = trajinList_.begin();
                                   traj != trajinList_.end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail 
    if ( (*traj)->BeginTraj(showProgress_) ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->TrajFilename().full());
      break;
    }
    // Set current parm from current traj.
    Topology* CurrentParm = (*traj)->TrajParm();
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != CurrentParm->Pindex());
#   ifdef TIMER
    setup_time.Start();
#   endif
    // If Parm has changed or trajectory velocity status has changed,
    // reset the frame.
    if (parmHasChanged || (TrajFrame.HasVelocity() != (*traj)->HasVelocity()))
      TrajFrame.SetupFrameV(CurrentParm->Atoms(), (*traj)->HasVelocity(), 
                            (*traj)->NreplicaDimension());

    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set active reference for this parm
      CurrentParm->SetReferenceCoords( refFrames_.ActiveReference() );
      // Set up actions for this parm
      if (actionList_.SetupActions( &CurrentParm )) {
        mprintf("WARNING: Could not set up actions for %s: skipping.\n",
                CurrentParm->c_str());
        continue;
      }
      lastPindex = CurrentParm->Pindex();
    }
#   ifdef TIMER
    setup_time.Stop();
#   endif
    // Loop over every Frame in trajectory
    (*traj)->PrintInfoLine();
#   ifdef TIMER
    trajin_time.Start();
    bool readMoreFrames = (*traj)->GetNextFrame(TrajFrame);
    trajin_time.Stop();
    while ( readMoreFrames )
#   else
    while ( (*traj)->GetNextFrame(TrajFrame) )
#   endif
    {
      // Since Frame can be modified by actions, save original and use CurrentFrame
      Frame* CurrentFrame = &TrajFrame;
      // Perform Actions on Frame
#     ifdef TIMER
      actions_time.Start();
#     endif
      bool suppress_output = actionList_.DoActions(&CurrentFrame, actionSet);
#     ifdef TIMER
      actions_time.Stop();
#     endif
      // Do Output
      if (!suppress_output) {
#       ifdef TIMER
        trajout_time.Start();
#       endif
        trajoutList_.Write(actionSet, CurrentParm, CurrentFrame);
#       ifdef TIMER
        trajout_time.Stop();
#       endif
      }
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
    readSets += (*traj)->NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  rprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);
# ifdef TIMER
  frames_time.Stop();
  mprintf("TIME: Trajectory processing occurred in %.4f seconds\n"
          "TIME: Avg. throughput= %.4f frames / second.\n"
          "TIME:\tTrajectory read took %.4f seconds (%.2f%% of processing).\n"
          "TIME:\tAction setup took %.4f seconds (%.2f%% of processing).\n"
          "TIME:\tAction frame processing took %.4f seconds (%.2f%% of processing).\n"
          "TIME:\tTrajectory output took %.4f seconds (%.2f%% of processing).\n",
          frames_time.Total(),
          (double)readSets / frames_time.Total(),
          trajin_time.Total(),  (trajin_time.Total()  / frames_time.Total() )*100.0, 
          setup_time.Total(),   (setup_time.Total()   / frames_time.Total() )*100.0, 
          actions_time.Total(), (actions_time.Total() / frames_time.Total() )*100.0,
          trajout_time.Total(), (trajout_time.Total() / frames_time.Total() )*100.0 );
# endif
  // Close output trajectories.
  trajoutList_.Close();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  mprintf("\nACTION OUTPUT:\n");
  actionList_.Print( );
# ifdef MPI
  // Sync DataSets across all threads. 
  //DSL_.SynchronizeData(); // NOTE: Disabled, trajs are not currently divided.
# endif
  // ========== A N A L Y S I S  P H A S E ==========
  mprintf("\nDATASETS:\n");
  if (!analysisList_.Empty()) {
    DSL_.List();
    RunAnalyses();
    // DEBUG: DataSets, post-Analysis
    mprintf("\nDATASETS AFTER ANALYSIS:\n");
  }
  DSL_.List();

  // ========== D A T A  W R I T E  P H A S E ==========
  // Print Datafile information
  DFL_.List();
  MasterDataFileWrite();

  return 0;
}

// CpptrajState::MasterDataFileWrite()
void CpptrajState::MasterDataFileWrite() {
  // Only Master does DataFile output
  if (worldrank==0) {
#   ifdef TIMER
    Timer datafile_time;
    datafile_time.Start();
#   endif
    DFL_.WriteAllDF();
#   ifdef TIMER
    datafile_time.Stop();
    mprintf("TIME: Data file write took %.4f seconds.\n", datafile_time.Total());
#   endif
  }
}

// CpptrajState::RunAnalyses()
int CpptrajState::RunAnalyses() {
  ++nrun_;
# ifdef TIMER
  Timer analysis_time;
  analysis_time.Start();
# endif
  int err = analysisList_.DoAnalyses();
# ifdef TIMER
  analysis_time.Stop();
  mprintf("TIME: Analyses took %.4f seconds.\n", analysis_time.Total());
# endif
  // If all Analyses completed successfully, clean up analyses.
  if ( err == 0) 
    analysisList_.Clear();
  return err;
}
