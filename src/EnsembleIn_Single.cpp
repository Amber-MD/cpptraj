#ifdef ENABLE_SINGLE_ENSEMBLE
#include "EnsembleIn_Single.h"
#include "TrajectoryFile.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
EnsembleIn_Single::EnsembleIn_Single() : eio_(0), ensembleSize_(0) {}

// DESTRUCTOR
EnsembleIn_Single::~EnsembleIn_Single() {
  if (eio_ != 0) {
    EndEnsemble();
    delete eio_;
  }
}

// EnsembleIn_Single::SetupEnsembleRead()
int EnsembleIn_Single::SetupEnsembleRead(FileName const& tnameIn, ArgList& argIn,
                                       Topology *tparmIn)
{
  if (eio_ != 0) delete eio_;
  // Set file name and topology pointer.
  if (SetTraj().SetNameAndParm(tnameIn, tparmIn)) return 1;
  // Detect file format
  TrajectoryFile::TrajFormatType tformat;
  if ( (eio_ = TrajectoryFile::DetectFormat( Traj().Filename(), tformat )) == 0 ) {
    mprinterr("Error: Could not determine trajectory %s format.\n", Traj().Filename().full());
    return 1;
  }
  eio_->SetDebug( debug_ );
  mprintf("\tReading '%s' as %s\n", Traj().Filename().full(), TrajectoryFile::FormatString(tformat));
  // Process ensemble args // TODO: Should be common to Ensemble?
  bool nosort = argIn.hasKey("nosort");
  // Process format-specific read args
  if (eio_->processReadArgs( argIn )) return 1;
  // Set up the format for reading and get the number of frames.
  int nframes = eio_->setupTrajin(Traj().Filename(), Traj().Parm());
  if (nframes == TrajectoryIO::TRAJIN_ERR) {
    mprinterr("Error: Could not set up %s for reading.\n", Traj().Filename().full());
    return 1;
  }
  if (debug_ > 0) {
    if (nframes != TrajectoryIO::TRAJIN_UNK)
      mprintf("\t'%s' contains %i frames.\n", Traj().Filename().base(), nframes);
    else
      mprintf("\t'%s' contains an unknown number of frames.\n",Traj().Filename().base());
  }
  // Set the start, stop, and offset args based on user input. Do some bounds
  // checking.
  if (SetTraj().Counter().CheckFrameArgs( nframes, argIn )) return 1;
  // Set trajectory coordinate info.
  cInfo_ = eio_->CoordInfo();
  // NOTE: ensembleSize_ is saved here as a shortcut. Should always equal whats in cInfo_
  // Determine if this trajectory actually contains an ensemble.
  // FIXME: Should check for < 2?
  ensembleSize_ = cInfo_.EnsembleSize();
  if (ensembleSize_ < 1) {
    mprinterr("Error: Cannot process single file ensemble with '%s'\n", 
              TrajectoryFile::FormatString(tformat));
    return 1;
  }
# ifdef MPI
  // Set up communicators
  if (Parallel::SetupComms( ensembleSize_ )) return 1;
# endif
  // If dimensions are present, assume search by indices, otherwise by temp.
  targetType_ = ReplicaInfo::NONE;
  if (cInfo_.ReplicaDimensions().Ndims() > 0)
    targetType_ = ReplicaInfo::INDICES;
  else if (cInfo_.HasTemp())
    targetType_ = ReplicaInfo::TEMP;
  else if (!nosort) {
    mprinterr("Error: Ensemble trajectory does not have indices or temperature.\n");
    return 1;
  }
  if (debug_ > 0)
    cInfo_.PrintCoordInfo( Traj().Filename().base(), Traj().Parm()->c_str() );
# ifdef MPI
  // This array will let each thread know who has what frame.
  frameidx_.resize( ensembleSize_ ); // TODO: Get rid of, should do all in TrajIO class.
# endif
  // Get a list of all temperatures/indices.
  TemperatureMap_.ClearMap();
  IndicesMap_.ClearMap();
  if (targetType_ == ReplicaInfo::TEMP || targetType_ == ReplicaInfo::INDICES )
  {
#   ifdef MPI
    FrameArray f_ensemble(1);
#   else
    FrameArray f_ensemble( ensembleSize_ );
#   endif
    f_ensemble.SetupFrames( Traj().Parm()->Atoms(), cInfo_ );
    if ( eio_->openTrajin() ) return 1;
    if ( eio_->readArray( Traj().Counter().Start(), f_ensemble ) ) return 1;
    eio_->closeTraj();
    if (targetType_ == ReplicaInfo::TEMP) {
      std::vector<double> allTemps( ensembleSize_, -1.0 );
#     ifdef MPI
      // Consolidate temperatures
      if (GatherTemperatures(f_ensemble[0].tAddress(), allTemps, EnsembleComm())) return 1;
#     else
      for (int en = 0; en != ensembleSize_; ++en)
        allTemps[en] = f_ensemble[en].Temperature();
#     endif
      if (SetTemperatureMap( allTemps )) return 1;
    } else if (targetType_ == ReplicaInfo::INDICES) {
      std::vector<RemdIdxType> allIndices( ensembleSize_ );
#     ifdef MPI
      // Consolidate replica indices
      if (GatherIndices(f_ensemble[0].iAddress(), allIndices, cInfo_.ReplicaDimensions().Ndims(),
                        EnsembleComm()))
        return 1;
#     else
      for (int en = 0; en != ensembleSize_; ++en)
        allIndices[en] = f_ensemble[en].RemdIndices();
#     endif
      if (SetIndicesMap( allIndices)) return 1;
    }
  }  

  return 0;
}

int EnsembleIn_Single::BeginEnsemble() {
  // Open the trajectory
  if (eio_->openTrajin()) {
    mprinterr("Error: Could not open %s\n",Traj().Filename().base());
    return 1;
  }
  // Initialize counter.
  SetTraj().Counter().Begin();
  return 0;
}

void EnsembleIn_Single::EndEnsemble() {
  if (eio_ != 0) {
    eio_->closeTraj();
#   ifdef MPI
#   ifdef TIMER
    total_mpi_allgather_ += mpi_allgather_timer_.Total();
    total_mpi_sendrecv_  += mpi_sendrecv_timer_.Total();
#   endif
#   endif
  }
}


int EnsembleIn_Single::ReadEnsemble(int currentFrame, FrameArray& f_ensemble,
                                  FramePtrArray& f_sorted )
{
  badEnsemble_ = false;
  // Read in all replicas.
  if ( eio_->readArray( currentFrame, f_ensemble ) ) return 1;
# ifdef MPI
  int ensembleFrameNum = 0;
  if (targetType_ != ReplicaInfo::NONE) {
    int my_idx;
    if (targetType_ == ReplicaInfo::TEMP)
      my_idx = TemperatureMap_.FindIndex( f_ensemble[0].Temperature() );
    else if (targetType_ == ReplicaInfo::INDICES)
      my_idx = IndicesMap_.FindIndex( f_ensemble[0].RemdIndices() );
#   ifdef TIMER
    mpi_allgather_timer_.Start();
#   endif
    // TODO: Put this in Traj_NcEnsemble
    if (EnsembleComm().AllGather( &my_idx, 1, MPI_INT, &frameidx_[0])) {
      rprinterr("Error: Gathering frame indices.\n");
      badEnsemble_ = true;
      return 0; // TODO: Better parallel error check
    }
    for (int i = 0; i != ensembleSize_; i++)
      if (frameidx_[i] == -1) {
        badEnsemble_ = true;
        break;
      }
#   ifdef TIMER
    mpi_allgather_timer_.Stop();
    mpi_sendrecv_timer_.Start();
#   endif
    if (!badEnsemble_) {
      for (int sendrank = 0; sendrank != ensembleSize_; sendrank++) {
        int recvrank = frameidx_[sendrank];
        if (sendrank != recvrank) {
          if (sendrank == Member())
            f_ensemble[0].SendFrame( recvrank, EnsembleComm() );
          else if (recvrank == Member()) {
            f_ensemble[1].RecvFrame( sendrank, EnsembleComm() );
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
# else
  int fidx;
  for (int i = 0; i != ensembleSize_; i++) {
    if (targetType_ == ReplicaInfo::TEMP)
      fidx = TemperatureMap_.FindIndex( f_ensemble[i].Temperature() );
    else if (targetType_ == ReplicaInfo::INDICES)
      fidx = IndicesMap_.FindIndex( f_ensemble[i].RemdIndices() );
    else // NONE
      fidx = i;
    if ( fidx == -1 )
      badEnsemble_ = true;
    else
      f_sorted[fidx] = &f_ensemble[i];
  }
# endif
  return 0;
}

void EnsembleIn_Single::EnsembleInfo(int showExtended) const {
  mprintf("'%s' (REMD ensemble size %i) ",Traj().Filename().base(), ensembleSize_); 
  eio_->Info();
  mprintf(", Parm %s",Traj().Parm()->c_str());
  if (cInfo_.HasBox()) mprintf(" (%s box)", cInfo_.TrajBox().TypeName());
  if (showExtended==1) Traj().Counter().PrintFrameInfo();
  if (debug_>0)
    mprintf(", %i atoms, Box %i",Traj().Parm()->Natom(),(int)cInfo_.HasBox());
}
#endif
