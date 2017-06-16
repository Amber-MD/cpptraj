#include <cstdio> // sscanf
#include "Traj_CharmmRestart.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

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
  if (infile_.IsOpen()) infile_.CloseFile();
}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
// TODO should read past everything up to coords or box
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
  fname_ = fname;
  if (infile_.OpenFileRead( fname_ )) return TRAJIN_ERR;
  // Read past first line. No error checks here, assume ID has taken care of that.
  infile_.Line();
  // Read past second blank line.
  infile_.Line();
  // Read number of title lines.
  const char* ptr = infile_.Line();
  int ntitle;
  sscanf(ptr, "%8i", &ntitle);
  mprintf("DEBUG: %i title lines.\n", ntitle);
  // Read title. Get rid of leading asterisks and trailing whitespace.
  std::string title;
  for (int i = 0; i != ntitle; i++) {
    ptr = infile_.Line();
    if (ptr == 0) {
      mprinterr("Error: Could not read expected number of title lines.\n");
      return TRAJIN_ERR;
    }
    int offset;
    if (ptr[0] != '*')
      offset = 0;
    else
      offset = 2;
    title.append( NoTrailingWhitespace(ptr+offset) + " " );
  }
  mprintf("DEBUG: TITLE:\n%s\n", title.c_str());
  SetTitle( title );

  // Seek down to next relevant section; !CRYSTAL or !NATOM
  Box cbox;
  while (ptr != 0 && ptr[0] != ' ' && ptr[1] != '!')
    ptr = infile_.Line();
  if (ptr[2] == 'C' && ptr[3] == 'R' && ptr[4] == 'Y') {
    // Has unit cell information. Read the shape matrix.
    // NOTE: Seems that the scientific exponent rep in Fortran
    //       can sometimes come out as 'D', which confuses
    //       sscanf, so replace that.
    // FIXME check for old version?
    ptr = infile_.Line();
    if (ptr == 0) return TRAJIN_ERR;
    double* bp = cbox.boxPtr();
    char buff[133];
    buff[132] = '\0';
    unsigned int idx = 0;
    for (const char* p = ptr; *p != '\0'; ++p, ++idx) {
      if (*p == 'D')
        buff[idx] = 'E';
      else
        buff[idx] = *p;
    }
    ptr = infile_.Line();
    if (ptr == 0) return TRAJIN_ERR;
    for (const char* p = ptr; *p != '\0'; ++p, ++idx) {
      if (*p == 'D')
        buff[idx] = 'E';
      else
        buff[idx] = *p;
    }
    sscanf(buff, "%22lE%22lE%22lE%22lE%22lE%22lE",
           bp, bp+1, bp+2, bp+3, bp+4, bp+5);
    mprintf("DEBUG: Shape Matrix: %g %g %g %g %g %g\n",
            bp[0], bp[1], bp[2], bp[3], bp[4], bp[5]);
  }

  closeTraj();
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
