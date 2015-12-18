#include "EnsembleIn_Multi.h"
#include "CpptrajStdio.h"
#include "DataFile.h" // TODO remove
#include "StringRoutines.h" // integerToString TODO remove
#ifdef MPI
#  include "Parallel.h"
#endif

// EnsembleIn_Multi::SetupEnsembleRead()
int EnsembleIn_Multi::SetupEnsembleRead(FileName const& tnameIn, ArgList& argIn, Topology *tparmIn)
{
  // Set file name and topology pointer.
  if (SetTraj().SetNameAndParm(tnameIn, tparmIn)) return 1;
  REMDtraj_.ClearIOarray();
  // Check for deprecated args
  if (argIn.hasKey("remdout")) {
    mprinterr("%s", TrajIOarray::DEPRECATED_remdout);
    return 1;
  }
  // Process any remlog keywords so they are not processed by SetupTrajIO
  std::string remlog_name = argIn.GetStringKey("remlog");
  double remlog_nstlim    = argIn.getKeyDouble("nstlim", 1.0);
  double remlog_ntwx      = argIn.getKeyDouble("ntwx",   1.0);
  bool no_sort = argIn.hasKey("nosort");
  // CRDIDXARG: Parse out 'crdidx <indices list>' now so it is not processed
  //            by SetupTrajIO.
  ArgList crdidxarg;
  if (argIn.Contains("crdidx"))
    crdidxarg.SetList( "crdidx " + argIn.GetStringKey("crdidx"), " " );
  // Set up replica file names.
  if (REMDtraj_.SetupReplicaFilenames( tnameIn, argIn )) return 1;
  // Set up TrajectoryIO classes for all file names.
  if (REMDtraj_.SetupIOarray(argIn, SetTraj().Counter(), cInfo_, Traj().Parm())) return 1;
# ifdef MPI
  // FIXME currently for ensemble processing need same # threads as replicas.
  if ((int)REMDtraj_.size() != Parallel::World().Size()) {
    mprinterr("Error: For ensemble processing # of threads (%i) must equal # replicas (%zu)\n",
              Parallel::World().Size(), REMDtraj_.size());
    return 1;
  }
# endif
  // Unless nosort was specified, figure out target type
  if (no_sort)
    targetType_ = ReplicaInfo::NONE;
  else {
    if ( !remlog_name.empty() ) {
      // Sort according to remlog data.
      DataFile remlogFile;
      DataSetList tempDSL;
      // CRDIDXARG: TODO: Come up with a way to do this that doesnt require ArgLists.
      if (remlogFile.ReadDataIn( remlog_name, crdidxarg, tempDSL ) ||
          tempDSL.empty())
      {
        mprinterr("Error: Could not read remlog data.\n");
        return 1;
      }
      if (remlogFile.Type() != DataFile::REMLOG) {
        mprinterr("Error: remlog: File was not of type remlog.\n");
        return 1;
      }
      if ( REMDtraj_.size() != tempDSL[0]->Size() ) {
        mprinterr("Error: ensemble size %zu does not match remlog ensemble size %zu\n",
                  REMDtraj_.size(), tempDSL[0]->Size());
        return 1;
      }
      remlogData_ = *((DataSet_RemLog*)tempDSL[0]); // FIXME: This feels clunky. Can we read direct?
      targetType_ = ReplicaInfo::CRDIDX;
      remdFrameFactor_ = remlog_ntwx / remlog_nstlim;
      mprintf("\t%g exchanges for every trajectory frame written.\n", remdFrameFactor_);
      if (remdFrameFactor_ > 1.0)
        remdFrameOffset_ = (int)remdFrameFactor_ - 1;
      else
        remdFrameOffset_ = 0;
      mprintf("\tTrajectory frame 1 corresponds to exchange %i\n", remdFrameOffset_ + 1);
      int expectedTrajFrames = (int)((double)Traj().Counter().TotalFrames() * remdFrameFactor_);
      if ( expectedTrajFrames != remlogData_.NumExchange() ) {
        mprinterr("Error: expected length of REMD ensemble %i does not match # exchanges in remlog %i.\n",
                  expectedTrajFrames, remlogData_.NumExchange());
        return 1;
      }
    } else {
      // If dimensions are present index by replica indices, otherwise index
      // by temperature. 
      if (cInfo_.ReplicaDimensions().Ndims() > 0)
        targetType_ = ReplicaInfo::INDICES;
      else
        targetType_ = ReplicaInfo::TEMP;
    }
  }
# ifdef MPI
  // This array will let each thread know who has what frame.
  frameidx_.resize( REMDtraj_.size() );
# endif
  // Get a list of all temperatures/indices.
  TemperatureMap_.ClearMap();
  IndicesMap_.ClearMap();
  if (targetType_ == ReplicaInfo::TEMP || targetType_ == ReplicaInfo::INDICES )
  {
    Frame frameIn;
    frameIn.SetupFrameV( Traj().Parm()->Atoms(), cInfo_ );
    std::vector<double> allTemps;
    std::vector<RemdIdxType> allIndices;
    if (targetType_ == ReplicaInfo::TEMP)
      allTemps.assign(REMDtraj_.size(), -1.0);
    else if (targetType_ == ReplicaInfo::INDICES)
      allIndices.resize( REMDtraj_.size() );
#   ifdef MPI
    int err = 0;
    if (REMDtraj_[Parallel::World().Rank()]->openTrajin()) {
      rprinterr("Error: Opening %s\n", REMDtraj_.f_name(Parallel::World().Rank()));
      err = 1;
    }
    if (err == 0) {
      if (REMDtraj_[Parallel::World().Rank()]->readFrame( Traj().Counter().Start(), frameIn )) {
        rprinterr("Error: Reading %s\n", REMDtraj_.f_name(Parallel::World().Rank()));
        err = 2;
      }
      REMDtraj_[Parallel::World().Rank()]->closeTraj();
    }
    int total_error;
    Parallel::World().AllReduce( &total_error, &err, 1, MPI_INT, MPI_SUM );
    if (total_error != 0) {
      mprinterr("Error: Cannot setup ensemble trajectories.\n");
      return 1;
    }
    if (targetType_ == ReplicaInfo::TEMP) {
      if (GatherTemperatures(frameIn.tAddress(), allTemps))
        return 1;
    } else if (targetType_ == ReplicaInfo::INDICES) {
      if (GatherIndices(frameIn.iAddress(), allIndices, cInfo_.ReplicaDimensions().Ndims()))
        return 1;
    }
#   else
    for (unsigned int member = 0; member != REMDtraj_.size(); member++) {
      if ( REMDtraj_[member]->openTrajin() ) return 1;
      if ( REMDtraj_[member]->readFrame( Traj().Counter().Start(), frameIn ) ) return 1;
      REMDtraj_[member]->closeTraj();
      if (targetType_ == ReplicaInfo::TEMP)
        allTemps[member] = frameIn.Temperature();
      else if (targetType_ == ReplicaInfo::INDICES)
        allIndices[member] = frameIn.RemdIndices();
    }
#   endif
    if (targetType_ == ReplicaInfo::TEMP) {
      if (SetTemperatureMap(allTemps)) return 1;
    } else if (targetType_ == ReplicaInfo::INDICES) {
      if (SetIndicesMap(allIndices)) return 1;
    }
  }  // Otherwise NONE, no sorting

  return 0;
}

// EnsembleIn_Multi::ReadEnsemble()
int EnsembleIn_Multi::ReadEnsemble(int currentFrame, FrameArray& f_ensemble, 
                                 FramePtrArray& f_sorted )
{
  int fidx = 0;
  badEnsemble_ = 0;
  // Read in all replicas
  //mprintf("DBG: Ensemble frame %i:",currentFrame+1); // DEBUG
# ifdef MPI
  int repIdx = Parallel::World().Rank(); // for targetType==CRDIDX
  unsigned int member = 0;
  // Read REMDtraj for this rank
  if ( REMDtraj_[Parallel::World().Rank()]->readFrame( currentFrame, f_ensemble[0]) )
    return 1;
# else
  int repIdx = 0; // for targetType==CRDIDX
  for (unsigned int member = 0; member != REMDtraj_.size(); ++member)
  {
    if ( REMDtraj_[member]->readFrame( currentFrame, f_ensemble[member]) )
      return 1;
# endif
    if (targetType_ == ReplicaInfo::TEMP)
      fidx = TemperatureMap_.FindIndex( f_ensemble[member].Temperature() );
    else if (targetType_ == ReplicaInfo::INDICES)
      fidx = IndicesMap_.FindIndex( f_ensemble[member].RemdIndices() );
    else if (targetType_ == ReplicaInfo::CRDIDX) {
      int currentRemExchange = (int)((double)currentFrame * remdFrameFactor_) + remdFrameOffset_;
      //mprintf("DEBUG:\tTrajFrame#=%i  RemdExch#=%i\n", currentFrame+1, currentRemExchange+1);
      fidx = remlogData_.RepFrame( currentRemExchange, repIdx++ ).CoordsIdx() - 1;
      //mprintf("DEBUG:\tFrame %i\tPosition %u is assigned index %i\n", currentFrame, member, fidx);
    }
#   ifndef MPI
    else // NONE
      fidx = member; // Not needed when MPI
    if (fidx == -1)
      badEnsemble_ = 1;
    else
      f_sorted[fidx] = &f_ensemble[member];
  }
#   else 
  // If calculated index is not worldrank, coords need to be sent to rank fidx.
  //rprintf("Index=%i\n", fidx); // DEBUG
  int ensembleFrameNum = 0;
  if (targetType_ != ReplicaInfo::NONE) {
    // Each rank needs to know where to send its coords, and where to receive coords from.
#   ifdef TIMER
    mpi_allgather_timer_.Start();
#   endif
    if (Parallel::World().AllGather( &fidx, 1, MPI_INT, &frameidx_[0])) {
      rprinterr("Error: Gathering frame indices.\n");
      return 0; // TODO: Better parallel error check
    }
    for (unsigned int i = 0; i != REMDtraj_.size(); i++)
      if (frameidx_[i] == -1) {
        badEnsemble_ = 1;
        break;
      }
#   ifdef TIMER
    mpi_allgather_timer_.Stop();
    mpi_sendrecv_timer_.Start();
#   endif
    // LOOP: one sendrecv at a time.
    if (badEnsemble_ == 0) {
      for (int sendrank = 0; sendrank != (int)REMDtraj_.size(); sendrank++) {
        int recvrank = frameidx_[sendrank];
        if (sendrank != recvrank) {
          if (sendrank == Parallel::World().Rank())
            f_ensemble[0].SendFrame( recvrank ); 
          else if (recvrank == Parallel::World().Rank()) {
            f_ensemble[1].RecvFrame( sendrank );
            // Since a frame was received, indicate position 1 should be used
            ensembleFrameNum = 1; 
          }
        }
        //else rprintf("SEND RANK == RECV RANK, NO COMM\n"); // DEBUG
      }
    }
#   ifdef TIMER
    mpi_sendrecv_timer_.Stop();
#   endif
  }
  f_sorted[0] = &f_ensemble[ensembleFrameNum];
  //rprintf("FRAME %i, FRAME RECEIVED= %i\n", currentFrame, ensembleFrameNum); // DEBUG 
#   endif
  return 0;
}

int EnsembleIn_Multi::BeginEnsemble() {
  mprintf("\tENSEMBLE: OPENING %zu REMD TRAJECTORIES\n", REMDtraj_.size());
# ifdef MPI
  // Open the trajectory this thread will be dealing with.
  if (REMDtraj_[Parallel::World().Rank()]->openTrajin()) {
    rprinterr("Error: Trajin_Multi::BeginTraj: Could not open replica %s\n",
              REMDtraj_.f_name(Parallel::World().Rank()));
    return 1;
  }
# else
  // Open the trajectories
  for (TrajIOarray::const_iterator rep = REMDtraj_.begin(); rep != REMDtraj_.end(); ++rep)
  {
    if ( (*rep)->openTrajin()) {
      mprinterr("Error: Could not open replica # %zu, '%s'\n",
                rep - REMDtraj_.begin(), REMDtraj_.f_name(rep-REMDtraj_.begin()));
      return 1;
    }
  }
# endif
  // Initialize counter.
  SetTraj().Counter().Begin();
  return 0;
}

void EnsembleIn_Multi::EndEnsemble() {
# ifdef MPI
  REMDtraj_[Parallel::World().Rank()]->closeTraj();
# ifdef TIMER
  total_mpi_allgather_ += mpi_allgather_timer_.Total();
  total_mpi_sendrecv_  += mpi_sendrecv_timer_.Total();
# endif
# else
  for (TrajIOarray::const_iterator rep = REMDtraj_.begin(); rep != REMDtraj_.end(); ++rep)
    (*rep)->closeTraj();
# endif
}

void EnsembleIn_Multi::EnsembleInfo(int showExtended) const {
  mprintf("Trajectory ensemble (%u total), lowest replica '%s'", REMDtraj_.size(),
          Traj().Filename().base());
  if (showExtended == 1) Traj().Counter().PrintFrameInfo();
  mprintf("\n");
  if (debug_ > 0) REMDtraj_.PrintIOinfo();
  if ( targetType_ == ReplicaInfo::INDICES )
    mprintf("\tProcessing ensemble using replica indices\n");
  else if ( targetType_ == ReplicaInfo::TEMP )
    mprintf("\tProcessing ensemble using replica temperatures\n");
  else if ( targetType_ == ReplicaInfo::CRDIDX )
    mprintf("\tProcessing ensemble using remlog data, sorting by coordinate index.\n");
  else // NONE 
    mprintf("\tNot sorting ensemble.\n");
  if (debug_ > 0) PrintReplicaInfo();
}

/** CRDIDXARG:
  * \return A string containing the coordinate indices (comma separated) of the
  *         final exchange in remlogData_.
  */
std::string EnsembleIn_Multi::FinalCrdIndices() const {
  if (remlogData_.Empty()) return std::string();
  std::string arg("crdidx ");
  int finalExchg = remlogData_.NumExchange() - 1;
  for (unsigned int rep = 0; rep < remlogData_.Size(); rep++) {
    if (rep > 0) arg += ",";
    arg += ( integerToString( remlogData_.RepFrame(finalExchg, rep).CoordsIdx() ) );
  }
  return arg;
}
