#include <cstdio> // sscanf
#include "Traj_CharmmRestart.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
Traj_CharmmRestart::Traj_CharmmRestart() {}

/** Identify trajectory format. File should be setup for READ */
bool Traj_CharmmRestart::ID_TrajFormat(CpptrajFile& fileIn) {
  if (fileIn.OpenFile()) return false;
  bool isCharmmRestart = false;
  const char* ptr = fileIn.NextLine();
  // Line begins: (A4,2I6,
  if (ptr != 0 && ptr[0] == 'R' && ptr[1] == 'E' && ptr[2] == 'S' && ptr[3] == 'T') {
    // Check if we can read 2 integers.
    int vernum, ldyna;
    if (sscanf(ptr+4, "%6i%6i", &vernum, &ldyna) == 2) {
      // Next line should be blank
      ptr = fileIn.NextLine();
      if (ptr != 0 && (ptr[0] == '\n' || ptr[0] == '\r')) {
        // Next line is I8 !NTITLE
        ptr = fileIn.NextLine();
        if (ptr != 0 && ptr[8] == ' ' && ptr[9] == '!' && ptr[10] == 'N') {
          isCharmmRestart = true;
          mprintf("DEBUG: Charmm restart file (%i): %s", vernum, ptr);
        }
      }
    }
  }
  fileIn.CloseFile();
  return isCharmmRestart;
}

/** Print trajectory info to stdout. */
void Traj_CharmmRestart::Info() {
  mprintf("is a CHARMM restart file");
}

/** Close file. */
void Traj_CharmmRestart::closeTraj() {

}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int Traj_CharmmRestart::openTrajin() {

  return 0;
}

/** Read help */
void Traj_CharmmRestart::ReadHelp() {

}

/** Process read arguments. */
int Traj_CharmmRestart::processReadArgs(ArgList& argIn) {

  return 0;
}

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int Traj_CharmmRestart::setupTrajin(FileName const& fname, Topology* trajParm)
{

  return TRAJIN_ERR;
}

/** Read specified trajectory frame. */
int Traj_CharmmRestart::readFrame(int set, Frame& frameIn) {

  return 0;
}

/** Read velocities from specified frame. */
int Traj_CharmmRestart::readVelocity(int set, Frame& frameIn) {

  return 0;
}

/** Read forces from specified frame. */
int Traj_CharmmRestart::readForce(int set, Frame& frameIn) {

  return 0;
}

// -----------------------------------------------------------------------------
/** Write help. */
void Traj_CharmmRestart::WriteHelp() {

}

/** Process write arguments. */
int Traj_CharmmRestart::processWriteArgs(ArgList& argIn) {

  return 0;
}

/** Set up trajectory for write. */
int Traj_CharmmRestart::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{

  return 1;
}

/** Write specified trajectory frame. */
int Traj_CharmmRestart::writeFrame(int set, Frame const& frameOut) {

  return 0;
}

// =============================================================================
#ifdef MPI
/** Open trajectory for reading in parallel. */
int Traj_CharmmRestart::parallelOpenTrajin(Parallel::Comm const& commIn) {
  return 1;
}

/** Open trajectory for writing in parallel. */
int Traj_CharmmRestart::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return 1;
}

/** Set up trajectory for write in parallel. */
int Traj_CharmmRestart::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{

  return 1;
}

/** Read frame in parallel. */
int Traj_CharmmRestart::parallelReadFrame(int set, Frame& frameIn) {

  return 1;
}

/** Write frame in parallel. */
int Traj_CharmmRestart::parallelWriteFrame(int set, Frame const& frameOut) {

  return 1;
}

/** Close trajectory in parallel. */
void Traj_CharmmRestart::parallelCloseTraj() {

}
#endif
