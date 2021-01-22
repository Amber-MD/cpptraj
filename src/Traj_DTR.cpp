#ifdef ENABLE_DTR
#include "Traj_DTR.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "Topology.h"
#include "vmdplugin/dtrplugin.hxx"

using namespace desres::molfile;
/// CONSTRUCTOR
Traj_DTR::Traj_DTR() :
  DTR_(0),
  fbuffer_(0),
  bufsize_(0)
{}

/// DESTRUCTOR
Traj_DTR::~Traj_DTR() {
  if (DTR_ != 0) delete DTR_;
  if (fbuffer_ != 0) delete[] fbuffer_;
}

/** Identify trajectory format. File should be setup for READ */
bool Traj_DTR::ID_TrajFormat(CpptrajFile& fileIn) {
  // For backwards compatibility with VMD
  if (fileIn.Filename().Base() == "clickme.dtr") return true;
  if (fileIn.OpenFile()) return false;
  unsigned char buffer[4];
  if (fileIn.Read(buffer, 4) != 4) return false;
  fileIn.CloseFile();
  if (buffer[0] != 'D' ||
      buffer[1] != 'E' ||
      buffer[2] != 'S' ||
      buffer[3] != 'M')
    return false;
  return true;
}

/** Print trajectory info to stdout. */
void Traj_DTR::Info() {
  mprintf("is a Desmond DTR trajectory");
}

/** Close file. */
void Traj_DTR::closeTraj() {

}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int Traj_DTR::openTrajin() {

  return 0;
}

/** Read help */
void Traj_DTR::ReadHelp() {

}

/** Process read arguments. */
int Traj_DTR::processReadArgs(ArgList& argIn) {

  return 0;
}

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int Traj_DTR::setupTrajin(FileName const& fname, Topology* trajParm)
{
  if (DTR_ != 0) delete DTR_;
  DTR_ = 0;
  if (fbuffer_ != 0) delete[] fbuffer_;
  fbuffer_ = 0;

  std::string initName;
  // check fot .stk file
  if (StkReader::recognizes(fname.full())) {
    DTR_ = new StkReader;
    initName = fname.Full();
  } else {
    DTR_ = new DtrReader;
    // DTR should be initialized with the dir name
    initName = fname.DirPrefix();
  }

  if (debug_ > 0) mprintf("DEBUG: initName= %s\n", initName.c_str());

  if (!DTR_->init( initName.c_str() )) {
    mprinterr("Error: DTR init failed.\n");
    delete DTR_;
    DTR_ = 0;
    return TRAJIN_ERR;
  }

  // Check number of atoms
  if ( (int)DTR_->natoms != trajParm->Natom() ) {
    mprinterr("Error: # of atoms in DTR (%u) != # atoms in associated topology (%i)\n",
              DTR_->natoms, trajParm->Natom());
    return TRAJIN_ERR;
  }

  ssize_t nframes = DTR_->size();

  if (debug_ > 0) mprintf("DEBUG: %zd frames.\n", nframes);
  if (nframes < 1) {
    mprinterr("Error: No frames detected in DTR trajectory.\n");
    return TRAJIN_ERR;
  }

  bool has_vel = DTR_->has_velocities();
  // Allocate float buffer
  bufsize_ = 3 * (size_t)trajParm->Natom();
  if (has_vel)
    bufsize_ *= 2;
  fbuffer_ = new float[ bufsize_ ];

  // Set Coordinate info.
  // Read the first frame to get the box info. A bit wasteful but gets
  // the job done for now.
  // NOTE: DTR seems to always have box? Always orthogonal?
  Box tmpBox;
  molfile_timestep_t Tstep;
  Tstep.coords = fbuffer_;
  if (has_vel)
    Tstep.velocities = fbuffer_ + bufsize_;
  int ret = DTR_->frame(0, &Tstep);
  // -1 is EOF
  // 0 is success
  if ( ret != 0 ) {
    mprinterr("Error: Could not read first frame of DTR during setup.\n");
    return 1;
  }

  tmpBox.SetupFromXyzAbg( Tstep.A, Tstep.B, Tstep.C,
                          Tstep.alpha, Tstep.beta, Tstep.gamma );

  SetCoordInfo( CoordinateInfo(tmpBox, has_vel, false, true) );
  
  return nframes;
}

static inline void FloatToDouble(double* darray, const float* farray, size_t asize)
{
  for (size_t idx = 0; idx != asize; idx++)
    darray[idx] = (double)farray[idx];
}

/** Read specified trajectory frame. */
int Traj_DTR::readFrame(int set, Frame& frameIn) {
  molfile_timestep_t Tstep;
  Tstep.coords = fbuffer_;
  if (CoordInfo().HasVel())
    Tstep.velocities = fbuffer_ + bufsize_;

  int ret = DTR_->frame(set, &Tstep);
  // -1 is EOF
  // 0 is success
  if ( ret != 0 ) return 1;
  // Convert float to double
  FloatToDouble(frameIn.xAddress(), Tstep.coords, bufsize_);
  if (CoordInfo().HasVel())
    FloatToDouble(frameIn.vAddress(), Tstep.velocities, bufsize_);

  // Set box
  frameIn.ModifyBox().AssignFromXyzAbg( Tstep.A, Tstep.B, Tstep.C,
                                        Tstep.alpha, Tstep.beta, Tstep.gamma );

  // Time
  frameIn.SetTime( Tstep.physical_time );

  return 0;
}

/** Read velocities from specified frame. */
int Traj_DTR::readVelocity(int set, Frame& frameIn) {

  return 0;
}

/** Read forces from specified frame. */
int Traj_DTR::readForce(int set, Frame& frameIn) {

  return 0;
}

// -----------------------------------------------------------------------------
/** Write help. */
void Traj_DTR::WriteHelp() {

}

/** Process write arguments. */
int Traj_DTR::processWriteArgs(ArgList& argIn, DataSetList const& DSLin) {

  return 0;
}

/** Set up trajectory for write. */
int Traj_DTR::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{

  return 1;
}

/** Write specified trajectory frame. */
int Traj_DTR::writeFrame(int set, Frame const& frameOut) {

  return 0;
}

// =============================================================================
#ifdef MPI
/** Open trajectory for reading in parallel. */
int Traj_DTR::parallelOpenTrajin(Parallel::Comm const& commIn) {
  return 1;
}

/** Open trajectory for writing in parallel. */
int Traj_DTR::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return 1;
}

/** Set up trajectory for write in parallel. */
int Traj_DTR::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{

  return 1;
}

/** Read frame in parallel. */
int Traj_DTR::parallelReadFrame(int set, Frame& frameIn) {

  return 1;
}

/** Write frame in parallel. */
int Traj_DTR::parallelWriteFrame(int set, Frame const& frameOut) {

  return 1;
}

/** Close trajectory in parallel. */
void Traj_DTR::parallelCloseTraj() {

}
#endif
#endif /* ENABLE_DTR */
