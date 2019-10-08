#ifndef NO_TNGFILE
#include <cmath>
#include "Traj_TNG.h"
#include "CpptrajStdio.h"
#include "FileName.h"
#include "Topology.h"

/// CONSTRUCTOR
Traj_TNG::Traj_TNG() {}

/** Identify trajectory format. File should be setup for READ */
bool Traj_TNG::ID_TrajFormat(CpptrajFile& fileIn) {

  return false;
}

/** Print trajectory info to stdout. */
void Traj_TNG::Info() {
  mprintf("is a <type>");
}

/** Close file. */
void Traj_TNG::closeTraj() {

}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int Traj_TNG::openTrajin() {

  return 0;
}

/** Read help */
void Traj_TNG::ReadHelp() {

}

/** Process read arguments. */
int Traj_TNG::processReadArgs(ArgList& argIn) {

  return 0;
}

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int Traj_TNG::setupTrajin(FileName const& fname, Topology* trajParm)
{
  tng_function_status stat = tng_util_trajectory_open(fname.full(), 'r', &traj_);
  if (stat != TNG_SUCCESS) {
    mprinterr("Error: Could not open TNG file '%s'\n", fname.full());
    return TRAJIN_ERR;
  }

  // Get number of atoms
  tng_num_particles_get(traj_, &tngatoms_);
  if (tngatoms_ != (int64_t)trajParm->Natom()) {
    mprinterr("Error: Number of atoms in TNG file (%li) does not match number\n"
              "Error:  of atoms in associated topology (%i)\n",
               tngatoms_, trajParm->Natom());
    return TRAJIN_ERR;
  }

  // Get number of frames
  int64_t nframes;
  tng_num_frames_get(traj_, &nframes);
  mprintf("\tTNG file has %li frames.\n", nframes);

  // Get the exponential distance scaling factor
  int64_t tngexp;
  if (tng_distance_unit_exponential_get(traj_, &tngexp) != TNG_SUCCESS) {
    mprinterr("Error: Could not get distance scaling exponential from TNG.\n");
    return TRAJIN_ERR;
  }
  mprintf("\tTNG exponential: %li\n", tngexp);
  switch (tngexp) {
    case 9  :
      // Input is in nm. Convert to Angstroms.
      //tngfac_ = 1.0;
      tngfac_ = 0.1;
      break;
    case 10 :
      // Input is in Angstroms.
      //tngfac_ = 10.0;
      tngfac_ = 1.0;
      break;
    default :
      // Convert to Angstroms.
      tngfac_ = pow(10.0, tngexp + 10); break;
  }
  mprintf("\tTNG distance scaling factor: %g\n", tngfac_);

  return (int)nframes;
}

/** Read specified trajectory frame. */
int Traj_TNG::readFrame(int set, Frame& frameIn) {

  return 0;
}

/** Read velocities from specified frame. */
int Traj_TNG::readVelocity(int set, Frame& frameIn) {

  return 0;
}

/** Read forces from specified frame. */
int Traj_TNG::readForce(int set, Frame& frameIn) {

  return 0;
}

// -----------------------------------------------------------------------------
/** Write help. */
void Traj_TNG::WriteHelp() {

}

/** Process write arguments. */
int Traj_TNG::processWriteArgs(ArgList& argIn, DataSetList const& DSLin) {

  return 0;
}

/** Set up trajectory for write. */
int Traj_TNG::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{

  return 1;
}

/** Write specified trajectory frame. */
int Traj_TNG::writeFrame(int set, Frame const& frameOut) {

  return 0;
}

// =============================================================================
#ifdef MPI
/** Open trajectory for reading in parallel. */
int Traj_TNG::parallelOpenTrajin(Parallel::Comm const& commIn) {
  return 1;
}

/** Open trajectory for writing in parallel. */
int Traj_TNG::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return 1;
}

/** Set up trajectory for write in parallel. */
int Traj_TNG::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{

  return 1;
}

/** Read frame in parallel. */
int Traj_TNG::parallelReadFrame(int set, Frame& frameIn) {

  return 1;
}

/** Write frame in parallel. */
int Traj_TNG::parallelWriteFrame(int set, Frame const& frameOut) {

  return 1;
}

/** Close trajectory in parallel. */
void Traj_TNG::parallelCloseTraj() {

}
#endif
#endif /* NO_TNGFILE */
