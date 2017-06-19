#include <cstdio> // sscanf
#include <cstdlib> // atof
#include "Traj_CharmmRestart.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "BufferedLine.h"

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
  // TODO Check read access?
}

/** Close file. */
void Traj_CharmmRestart::closeTraj() {
  if (inframe_.IsOpen()) inframe_.CloseFile();
}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int Traj_CharmmRestart::openTrajin() {
  if (inframe_.OpenFile()) return 1;
  return 0;
}

/** Read help */
void Traj_CharmmRestart::ReadHelp() {

}

/** Process read arguments. */
int Traj_CharmmRestart::processReadArgs(ArgList& argIn) {

  return 0;
}

static inline int ErrEOF(int line) {
  mprinterr("Error: Unexpected end of file, line %i\n");
  return TrajectoryIO::TRAJIN_ERR;
}

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int Traj_CharmmRestart::setupTrajin(FileName const& fname, Topology* trajParm)
{
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) return TRAJIN_ERR;
  // Read past first line. No error checks here, assume ID has taken care of that.
  infile.Line();
  // Read past second blank line.
  infile.Line();
  // Read number of title lines.
  const char* ptr = infile.Line();
  int ntitle;
  sscanf(ptr, "%8i", &ntitle);
  mprintf("DEBUG: %i title lines.\n", ntitle);
  // Read title. Get rid of leading asterisks and trailing whitespace.
  std::string title;
  for (int i = 0; i != ntitle; i++) {
    ptr = infile.Line();
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
  cbox_.SetNoBox();
  while (ptr != 0 && ptr[1] != '!')
    ptr = infile.Line();
  if (ptr == 0) return ErrEOF(infile.LineNumber());
  if (ptr[2] == 'C' && ptr[3] == 'R' && ptr[4] == 'Y') {
    // Has unit cell information. Read the shape matrix.
    // NOTE: Seems that the scientific exponent rep in Fortran
    //       can sometimes come out as 'D', which confuses
    //       sscanf, so replace that.
    // FIXME check for old version?
    ptr = infile.Line();
    if (ptr == 0) return ErrEOF(infile.LineNumber());
    double bs[6];
    char buff[133];
    buff[132] = '\0';
    unsigned int idx = 0;
    for (const char* p = ptr; *p != '\0'; ++p, ++idx) {
      if (*p == 'D')
        buff[idx] = 'E';
      else
        buff[idx] = *p;
    }
    ptr = infile.Line();
    if (ptr == 0) return ErrEOF(infile.LineNumber());
    for (const char* p = ptr; *p != '\0'; ++p, ++idx) {
      if (*p == 'D')
        buff[idx] = 'E';
      else
        buff[idx] = *p;
    }
    sscanf(buff, "%22lE%22lE%22lE%22lE%22lE%22lE",
           bs, bs+1, bs+2, bs+3, bs+4, bs+5);
    mprintf("DEBUG: Shape Matrix: %g %g %g %g %g %g\n",
            bs[0], bs[1], bs[2], bs[3], bs[4], bs[5]);
    double bp[6];
    Box::ShapeToUcell( bp, bs );
    cbox_.SetBox( bp );
    mprintf("DEBUG: Unit cell: %g %g %g %g %g %g\n",
            cbox_.BoxX(), cbox_.BoxY(), cbox_.BoxZ(),
            cbox_.Alpha(), cbox_.Beta(), cbox_.Gamma());
    // Seek down to !NATOM
    while (ptr != 0 && ptr[1] != '!')
      ptr = infile.Line();
    if ( ptr == 0 ) return ErrEOF(infile.LineNumber());
    if (ptr[2] != 'N' || ptr[3] != 'A' || ptr[4] != 'T') {
      mprinterr("Error: !NATOM section not found.\n");
      return TRAJIN_ERR;
    }
  }

  ptr = infile.Line();
  int natom = 0;
  sscanf(ptr, "%12i", &natom);
  mprintf("DEBUG: %i atoms.\n", natom);
  if (natom != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in restart (%i) does not match # in\n", natom);
    mprinterr("Error:  associated topology %s (%i)\n", trajParm->c_str(),
              trajParm->Natom());
  }
  ncoord_ = natom * 3;
  // TODO temperature, time
  infile.CloseFile();
  // Box, vel, temp, time
  SetCoordInfo( CoordinateInfo(cbox_, true, false, false) );
  // Setup input frame
  if (inframe_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  inframe_.SetupFrameBuffer(ncoord_, 22, 3);
  
  // Only 1 frame
  return 1;
}

/** Read a block of coordinates laid out XYZ0 XYZ1 ... */
int Traj_CharmmRestart::ReadXYZ(double* xAddress) {
  if (inframe_.ReadFrame()) return 1;
  inframe_.BufferBegin();
  double* xptr = xAddress;
  char buf[23];
  buf[22] = '\0';
  for (int elt = 0; elt != ncoord_; elt++, xptr++) {
    const char* ptr = inframe_.NextElement();
    std::copy(ptr, ptr+22, buf);
    buf[18] = 'E';
    *xptr = atof( buf );
  }

  return 0;
}

/** Read specified trajectory frame. */
int Traj_CharmmRestart::readFrame(int set, Frame& frameIn) {
  // Seek down to ' !XOLD'
  const char* ptr = inframe_.NextLine();
  while (ptr != 0 && (ptr[0] != ' ' || ptr[1] != '!' ||
                      ptr[2] != 'X' || ptr[3] != 'O'))
    ptr = inframe_.NextLine();
  mprintf("DEBUG: %s\n", ptr);
  // Read coords 
  ReadXYZ(frameIn.xAddress());
  // Read velocities
  if (readVelocity(set, frameIn)) return 1;
  // Set box from what was read in setupTrajin().
  frameIn.SetBox( cbox_ );
  return 0;
}

/** Read velocities from specified frame. */
int Traj_CharmmRestart::readVelocity(int set, Frame& frameIn) {
  // Seek down to ' !VX'
  const char* ptr = inframe_.NextLine();
  while (ptr != 0 && (ptr[0] != ' ' || ptr[1] != '!' ||
                      ptr[2] != 'V' || ptr[3] != 'X'))
    ptr = inframe_.NextLine();
  mprintf("DEBUG: %s\n", ptr);
  // Read Velocities 
  ReadXYZ(frameIn.vAddress());
  return 0;
}

/** Read forces from specified frame. */
int Traj_CharmmRestart::readForce(int set, Frame& frameIn) {

  return 1;
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
