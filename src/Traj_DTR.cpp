#include "Traj_DTR.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
Traj_DTR::Traj_DTR() :
  DTR_(0)
{}

/** Identify trajectory format. File should be setup for READ */
bool Traj_DTR::ID_TrajFormat(CpptrajFile& fileIn) {

  return false;
}

/** Print trajectory info to stdout. */
void Traj_DTR::Info() {
  mprintf("is a <type>");
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

  return TRAJIN_ERR;
}

/** Read specified trajectory frame. */
int Traj_DTR::readFrame(int set, Frame& frameIn) {

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
