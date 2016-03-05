#include <cstdio> // sscanf
#include "Traj_CharmmCor.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"

bool Traj_CharmmCor::ID_TrajFormat(CpptrajFile& fileIn) {
  // File must already be set up for read.
  if (fileIn.OpenFile()) return false;
  bool isCor = false;
  const char* ptr = fileIn.NextLine();
  // Must be at least 1 title line denoted with '*'
  if (ptr != 0 && *ptr == '*') {
    // Scan past all title lines
    while (ptr != 0 && *ptr == '*') ptr = fileIn.NextLine();
    if (ptr != 0) {
      // Next line must be # atoms ONLY
      int ibuf[2];
      if (sscanf(ptr, "%i %i", ibuf, ibuf+1) == 1)
        // make sure it was a valid integer
        isCor = (ibuf[0] > 0);
    }
  }
  fileIn.CloseFile();
  return isCor;
}

int Traj_CharmmCor::setupTrajin(FileName const& fname, Topology* trajParm)
{
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  if (file_.OpenFile()) return TRAJIN_ERR;
  // Read past header.
  const char* buffer = file_.NextLine();
  if (buffer == 0) return TRAJIN_ERR;
  // Use first title line as traj title.
  const char* ptr = buffer;
  while (*ptr != '\0' && (*ptr == '*' || *ptr == ' ')) ++ptr;
  SetTitle( NoTrailingWhitespace(std::string(ptr)) );
  // Advance past header
  while (buffer != 0 && buffer[0] == '*') buffer = file_.NextLine();
  // Should now be positioned at # atoms. Check for EXT keyword.
  ArgList atomLine(buffer);
  extendedFmt_ = atomLine.hasKey("EXT"); 
  corAtom_ = atomLine.getNextInteger(-1);
  mprintf("\tCOR file: %i atoms\n", corAtom_);
  if (corAtom_ < 1) {
    mprinterr("Error: No atoms in CHARMM COR file.\n");
    return TRAJIN_ERR;
  }
  if (corAtom_ > 99999)
    extendedFmt_ = true;
  if (corAtom_ != trajParm->Natom()) {
    mprinterr("Error: COR file has %i atoms, associated topology '%s' has %i\n",
              corAtom_, trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  file_.CloseFile();
  // Just 1 frame.
  return 1;
}

int Traj_CharmmCor::openTrajin() { 
  if (file_.OpenFile()) return 1;
  // Advance past header
  const char* buffer = file_.NextLine();
  while (buffer != 0 && buffer[0] == '*') buffer = file_.NextLine();
  if (buffer == 0) return 1;
  return 0;
}

void Traj_CharmmCor::closeTraj() { file_.CloseFile(); }

int Traj_CharmmCor::readFrame(int set, Frame& frameIn) {
  // Atom #, Res #,   Res name, At name, X, Y, Z,       segid,   res id, weight
  // Original:
  // I5,     I5, 1X,  A4, 1X,   A4,      3(F10.5),  1X, A4, 1X,  A4,     F10.5
  // Extended:
  // I10,    I10, 1X, A9, 1X,  A9,     3(F20.10), 1X, A9, 1X, A9,    F20.10
  // Should be positioned at first atom line.
  double* xptr = frameIn.xAddress();
  for (int at = 0; at != corAtom_; at++, xptr += 3) {
    const char* buffer = file_.NextLine();
    if (buffer == 0) {
      mprinterr("Error: Reading COR atom %i\n", at+1);
      return 1;
    }
    int ncrd;
    if (extendedFmt_)
      ncrd = sscanf(buffer, "%*5i%*5i%*5s%*5s%10lf%10lf%10lf\n", xptr, xptr+1, xptr+2);
    else
      ncrd = sscanf(buffer, "%*10i%*10i%*10s%*10s%20lf%20lf%20lf\n", xptr, xptr+1, xptr+2);
    if (ncrd != 3) {
      mprinterr("Error: Reading coordinates for COR atom %i (got %i)\n", at+1, ncrd);
      return 1;
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
int Traj_CharmmCor::setupTrajout(FileName const& fname, Topology* trajParm,
                               CoordinateInfo const& cInfoIn,
                               int NframesToWrite, bool append)
{
  return 1;
}

int Traj_CharmmCor::writeFrame(int set, Frame const& frameOut) { return 1; }

// -----------------------------------------------------------------------------
void Traj_CharmmCor::Info() {
  mprintf("is a CHARMM COR file");
  if (extendedFmt_) mprintf(" (extended format)");
}
