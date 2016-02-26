#include "Trajin_Single.h"
#include "TrajectoryFile.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Trajin_Single::Trajin_Single() : trajio_(0), velio_(0), frcio_(0) {}

// DESTRUCTOR
Trajin_Single::~Trajin_Single() {
  if (trajio_!=0) {
    EndTraj();
    delete trajio_;
  }
  if (velio_!=0) delete velio_;
}

// TODO: Should this take a FileName instead of string?
int Trajin_Single::SetupTrajRead(FileName const& tnameIn, ArgList& argIn, 
                                 Topology* tparmIn)
{
  if (trajio_ != 0) delete trajio_;
  if (velio_ != 0) delete velio_;
  if (frcio_ != 0) delete frcio_;
  // Set file name and topology pointer.
  if (SetTraj().SetNameAndParm(tnameIn, tparmIn)) return 1;
  // Detect file format
  TrajectoryFile::TrajFormatType tformat;
  if ( (trajio_ = TrajectoryFile::DetectFormat( Traj().Filename(), tformat )) == 0 ) {
    mprinterr("Error: Could not determine trajectory %s format.\n", Traj().Filename().full());
    return 1;
  }
  trajio_->SetDebug( debug_ );
  mprintf("\tReading '%s' as %s\n", Traj().Filename().full(), TrajectoryFile::FormatString(tformat));
  // Process format-specific read args
  if (trajio_->processReadArgs( argIn )) return 1;
  // Set up the format for reading and get the number of frames.
  int nframes = trajio_->setupTrajin(Traj().Filename(), Traj().Parm());
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
  cInfo_ = trajio_->CoordInfo();
  // Check if a separate mdvel file will be read
  if (argIn.Contains("mdvel")) {
    FileName mdvel_fname( argIn.GetStringKey("mdvel") );
    if (mdvel_fname.empty()) {
      mprinterr("Error: mdvel: Usage 'mdvel <velocity filename>'\n");
      return 1;
    }
    if ( !File::Exists( mdvel_fname ) ) {
      File::ErrorMsg( mdvel_fname.full() );
      return 1;
    }
    // Detect mdvel format
    if ( (velio_ = TrajectoryFile::DetectFormat( mdvel_fname, tformat )) == 0 ) {
      mprinterr("Error: Could not set up velocity file %s for reading.\n",mdvel_fname.full());
      return 1;
    }
    velio_->SetDebug( debug_ );
    // Set up the format for reading mdvel, get # of mdvel frames
    int vel_frames = velio_->setupTrajin(mdvel_fname.Full(), Traj().Parm());
    if (vel_frames != Traj().Counter().TotalFrames()) {
      mprinterr("Error: velocity file %s frames (%i) != traj file frames (%i)\n",
                mdvel_fname.full(), vel_frames, Traj().Counter().TotalFrames());
      return 1;
    }
    cInfo_.SetVelocity( true );
  }

  // TODO add in support for separate mdfrc file
  if (debug_ > 0)
    cInfo_.PrintCoordInfo( Traj().Filename().base(), Traj().Parm()->c_str() );
  return 0;
}

int Trajin_Single::BeginTraj() {
  // Open the trajectory
  if (trajio_->openTrajin()) {
    mprinterr("Error: Trajin_Single::BeginTraj: Could not open %s\n",Traj().Filename().base());
    return 1;
  }
  // Open mdvel file if present
  if (velio_!=0 && velio_->openTrajin()) {
    mprinterr("Error: Could not open mdvel file.\n");
    return 1;
  }
  // TODO open mdfrc file if present
  // Initialize counter.
  SetTraj().Counter().Begin();
  return 0;
}

void Trajin_Single::EndTraj() {
  trajio_->closeTraj();
  if (velio_ != 0) velio_->closeTraj();
  // TODO close mdfrc file if present
}

// Trajin_Single::ReadTrajFrame()
int Trajin_Single::ReadTrajFrame( int idx, Frame& frameIn ) {
  if (trajio_->readFrame(idx, frameIn))
    return 1;
  if (velio_ != 0 && velio_->readVelocity(idx, frameIn))
    return 1;
  if (frcio_ != 0 && frcio_->readForce(idx, frameIn))
    return 1;
  //printf("DEBUG:\t%s:  current=%i  target=%i\n",trajName,idx,targetSet);
  return 0;
}

// Trajin_Single::PrintInfo()
void Trajin_Single::PrintInfo(int showExtended) const {
  mprintf("'%s' ",Traj().Filename().base());
  trajio_->Info();
  mprintf(", Parm %s",Traj().Parm()->c_str());
  if (trajio_->CoordInfo().HasBox())
    mprintf(" (%s box)", trajio_->CoordInfo().TrajBox().TypeName());
  if (showExtended==1) Traj().Counter().PrintFrameInfo(); 
  if (debug_>0)
    mprintf(", %i atoms, Box %i",Traj().Parm()->Natom(),(int)trajio_->CoordInfo().HasBox());
  mprintf("\n");
  if (velio_!=0) {
    mprintf("\tMDVEL: ");
    velio_->Info();
    mprintf("\n");
  }
  if (frcio_!=0) {
    mprintf("\tMDFRC: ");
    frcio_->Info();
    mprintf("\n");
  }
}

#ifdef MPI
int Trajin_Single::ParallelBeginTraj( Parallel::Comm const& commIn ) {
  // Open the trajectory
  if (trajio_->parallelOpenTrajin( commIn )) {
    rprinterr("Error: Could not open '%s' in parallel.\n", Traj().Filename().base());
    return 1;
  }
  // Open mdvel file if present // TODO enable?
  if (velio_ != 0) {
    mprinterr("Error: separate velocity file not yet supported in parallel.\n");
    return 1;
  }
  //if (velio_!=0 && velio_->openTrajin()) {
  //  mprinterr("Error: Could not open mdvel file.\n");
  //  return 1;
  //}
  // TODO open mdfrc file if present
  // Initialize counter. // FIXME needed in parallel?
  SetTraj().Counter().Begin();
  return 0;
}

int Trajin_Single::ParallelReadTrajFrame(int idx, Frame& frameIn) {
  if (trajio_->parallelReadFrame(idx, frameIn))
    return 1;
  //if (velio_ != 0 && velio_->readVelocity(idx, frameIn))
  //  return 1;
  //if (frcio_ != 0 && frcio_->readForce(idx, frameIn))
  //  return 1;
  //printf("DEBUG:\t%s:  current=%i  target=%i\n",trajName,idx,targetSet);
  return 0;
}

void Trajin_Single::ParallelEndTraj() {
  trajio_->parallelCloseTraj();
  //if (velio_ != 0) velio_->closeTraj();
  // TODO close mdfrc file if present
}
#endif
