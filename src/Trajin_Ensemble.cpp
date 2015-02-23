#ifdef ENABLE_SINGLE_ENSEMBLE
#include "Trajin_Ensemble.h"
#include "CpptrajStdio.h"
#ifdef MPI
# include "MpiRoutines.h"
#  ifdef TIMER
double Trajin_Ensemble::total_mpi_allgather_ = 0.0;
double Trajin_Ensemble::total_mpi_sendrecv_ = 0.0;
#  endif
#endif

// CONSTRUCTOR
Trajin_Ensemble::Trajin_Ensemble() :
  targetType_(ReplicaInfo::NONE),
  eio_(0),
  ensembleSize_(0),
  trajIsOpen_(false),
  badEnsemble_(false)
{}

// DESTRUCTOR
Trajin_Ensemble::~Trajin_Ensemble() {
  EndTraj();
  if (eio_ != 0) delete eio_;
}

int Trajin_Ensemble::SetupTrajRead(std::string const& tnameIn, ArgList& argIn, Topology *tparmIn)
{
  // Require a filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: Trajin_Single: No filename given.\n");
    return 1;
  }
  // Check and set associated parm file
  if ( SetTrajParm( tparmIn ) ) return 1;
  // Check that file can be opened. 
  //if (!fileExists(tnameIn)) return 1;
  // Detect file format
  TrajFormatType tformat;
  if ( (eio_ = DetectFormat( tnameIn, tformat )) == 0 ) {
    mprinterr("Error: Could not determine trajectory %s format.\n", tnameIn.c_str());
    return 1;
  }
  eio_->SetDebug( debug_ );
  // Set trajectory filename
  SetTrajFileName( tnameIn, true );
  mprintf("\tReading '%s' as %s\n", TrajFilename().full(), TrajectoryFile::FormatString(tformat));
  bool nosort = argIn.hasKey("nosort");
  // Process format-specific read args
  if (eio_->processReadArgs( argIn )) return 1;
  // Set up the format for reading and get the number of frames.
  if (SetupTrajIO( tnameIn, *eio_, argIn )) return 1;
  // Set trajectory coordinate info.
  cInfo_ = eio_->CoordInfo();
  // NOTE: ensembleSize_ is saved here as a shortcut. Should always equal whats in cInfo_
  // Determine if this trajectory actually contains an ensemble.
  // FIXME: Should check for < 2?
  ensembleSize_ = cInfo_.EnsembleSize();
  if (ensembleSize_ < 1) {
    mprinterr("Error: Cannot process single file ensemble with '%s'\n", FormatString(tformat));
    return 1;
  }
  // Check how many frames will actually be read
  if (setupFrameInfo() == 0) return 1;
  // Check traj box info against parm box info
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
    Frame::PrintCoordInfo( TrajFilename().base(), TrajParm()->c_str(), cInfo_ );
  // FIXME: Should this ever be done here?
  TrajParm()->SetParmCoordInfo( cInfo_ );
  return 0;
}

int Trajin_Ensemble::BeginTraj(bool showProgress) {
  // Open the trajectory
  if (eio_->openTrajin()) {
    mprinterr("Error: Could not open %s\n",TrajFilename().base());
    return 1;
  }
  // Set progress bar, start and offset.
  PrepareForRead( showProgress );
  trajIsOpen_ = true;
  return 0;
}

void Trajin_Ensemble::EndTraj() {
  if (trajIsOpen_) {
    eio_->closeTraj();
    trajIsOpen_ = false;
  }
# ifdef MPI
# ifdef TIMER
  total_mpi_allgather_ += mpi_allgather_timer_.Total();
  total_mpi_sendrecv_  += mpi_sendrecv_timer_.Total();
# endif
# endif
}

#ifdef MPI
#ifdef TIMER
void Trajin_Ensemble::TimingData(double trajin_time) {
  if (total_mpi_allgather_ > 0.0 || total_mpi_sendrecv_ > 0.0) {
    double other_time = trajin_time - total_mpi_allgather_ - total_mpi_sendrecv_;
    rprintf("MPI_TIME:\tallgather: %.4f s (%.2f%%), sendrecv: %.4f s (%.2f%%), Other:  %.4f s (%.2f%%)\n",
            total_mpi_allgather_, (total_mpi_allgather_ / trajin_time)*100.0,
            total_mpi_sendrecv_,  (total_mpi_sendrecv_  / trajin_time)*100.0,
            other_time, (other_time / trajin_time)*100.0 );
  }
}
#endif
#endif

void Trajin_Ensemble::PrintInfo(int showExtended) const {
  mprintf("'%s' (REMD ensemble size %i) ",TrajFilename().base(), ensembleSize_); 
  eio_->Info();
  mprintf(", Parm %s",TrajParm()->c_str());
  if (cInfo_.HasBox()) mprintf(" (%s box)", cInfo_.TrajBox().TypeName());
  if (showExtended==1) PrintFrameInfo();
  if (debug_>0)
    mprintf(", %i atoms, Box %i",TrajParm()->Natom(),(int)cInfo_.HasBox());
}

// -----------------------------------------------------------------------------
// Trajin_Ensemble::EnsembleInfo()
void Trajin_Ensemble::EnsembleInfo() const {
  if (targetType_ == ReplicaInfo::TEMP)
    PrintReplicaTmap( TemperatureMap_ );
  else if (targetType_ == ReplicaInfo::INDICES)
    PrintReplicaImap( IndicesMap_ );
}

// Trajin_Ensemble::EnsembleSetup()
int Trajin_Ensemble::EnsembleSetup( FrameArray& f_ensemble, FramePtrArray& f_sorted ) {
  // Allocate space to hold position of each incoming frame in replica space.
# ifdef MPI
  // Only two frames needed; one for reading, one for receiving.
  f_sorted.resize( 2 );
  f_ensemble.resize( 2 );
  // This array will let each thread know who has what frame.
  frameidx_.resize( ensembleSize_ ); // TODO: Get rid of, should do all in TrajIO class.
# else
  f_sorted.resize( ensembleSize_ );
  f_ensemble.resize( ensembleSize_ );
# endif 
  f_ensemble.SetupFrames( TrajParm()->Atoms(), cInfo_ );
  // Get a list of all temperatures/indices.
  TemperatureMap_.ClearMap();
  IndicesMap_.ClearMap();
  if (targetType_ == ReplicaInfo::TEMP || targetType_ == ReplicaInfo::INDICES )
  {
    if ( eio_->openTrajin() ) return 1;
    if ( eio_->readArray( Start(), f_ensemble ) ) return 1;
    eio_->closeTraj();
    if (targetType_ == ReplicaInfo::TEMP) {
      TemperatureMap_ = SetReplicaTmap(ensembleSize_, f_ensemble);
      if (TemperatureMap_.empty()) return 1;
    } else if (targetType_ == ReplicaInfo::INDICES) {
      IndicesMap_ = SetReplicaImap(ensembleSize_, cInfo_.ReplicaDimensions().Ndims(), f_ensemble);
      if (IndicesMap_.empty()) return 1;
    }
  }

  return 0;
}

// Trajin_Ensemble::ReadEnsemble()
int Trajin_Ensemble::ReadEnsemble(int currentFrame, FrameArray& f_ensemble,
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
    if (parallel_allgather( &my_idx, 1, PARA_INT, &frameidx_[0], 1, PARA_INT)) {
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
          if (sendrank == worldrank)
            f_ensemble[0].SendFrame( recvrank );
          else if (recvrank == worldrank) {
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
#endif
