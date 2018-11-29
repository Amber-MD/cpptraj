#include "Traj_XYZ.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

/// CONSTRUCTOR
Traj_XYZ::Traj_XYZ() :
  titleType_(NO_TITLE),
  ftype_(UNKNOWN),
  set_(0),
  fmt_(0)
{}

/** \param line1 First line. Will be set to title line if title is present, cleared otherwise.
  * \param line2 Second line.
  */
Traj_XYZ::Type Traj_XYZ::DetermineFormat(std::string& line1,
                                         std::string const& line2)
const
{
  std::string line = line1;
  // This line will either be a title line, atom XYZ, or XYZ
  RemoveAllWhitespace( line );
  if (!line.empty() && line[0] == '#') {
    line1 = line;
    // Ignore the title. Get the next line, which must be atom XYZ or XYZ
    line = line2;
  } else
    line1.clear();
  static const unsigned int tkSize = 64;
  char tk0[tkSize];
  char tk1[tkSize];
  char tk2[tkSize];
  char tk3[tkSize];
  int nscan = sscanf(line.c_str(), "%s %s %s %s", tk0, tk1, tk2, tk3);
  if (nscan == 4) {
    if ( validInteger( std::string(tk0) ) &&
         validDouble( std::string(tk1) ) &&
         validDouble( std::string(tk2) ) &&
         validDouble( std::string(tk3) ) )
      return ATOM_XYZ;
  } else if (nscan == 3) {
    if ( validDouble( std::string(tk0) ) &&
         validDouble( std::string(tk1) ) &&
         validDouble( std::string(tk2) ) )
      return XYZ;
  }
  return UNKNOWN;
}

/** Identify trajectory format. File should be setup for READ */
bool Traj_XYZ::ID_TrajFormat(CpptrajFile& fileIn) {
  // File must already be set up for read.
  if (fileIn.OpenFile()) return false;
  std::string line1 = fileIn.GetLine();
  std::string line2 = fileIn.GetLine();
  if ( DetermineFormat(line1, line2) == UNKNOWN ) return false;

  return true;
}

/** Print trajectory info to stdout. */
void Traj_XYZ::Info() {
  mprintf("is an XYZ trajectory");
}

/** Close file. */
void Traj_XYZ::closeTraj() {
  outfile_.CloseFile();
  infile_.CloseFile();
}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int Traj_XYZ::openTrajin() {
  if (infile_.OpenFile()) return 1;
  set_ = 0;
  return 0;
}

/** Read help */
void Traj_XYZ::ReadHelp() {

}

/** Process read arguments. */
int Traj_XYZ::processReadArgs(ArgList& argIn) {

  return 0;
}

const char* Traj_XYZ::FMT_XYZ_ = "%lf %lf %lf";

const char* Traj_XYZ::FMT_ATOM_XYZ_ = "%*i %lf %lf %lf";

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int Traj_XYZ::setupTrajin(FileName const& fname, Topology* trajParm)
{
  if (infile_.OpenFileRead( fname )) return TRAJIN_ERR;
  // Initial format determiniation.
  std::string line1 = infile_.GetLine();
  std::string line2 = infile_.GetLine();
  ftype_ = DetermineFormat(line1, line2);
  switch (ftype_) {
    case XYZ      : fmt_ = FMT_XYZ_; break;
    case ATOM_XYZ : fmt_ = FMT_ATOM_XYZ_; break;
    case UNKNOWN  :
      mprinterr("Internal Error: '%s' does not appear to be XYZ format anymore.\n",
                infile_.Filename().full());
      return TRAJIN_ERR;
  }
  if (line1.empty())
    titleType_ = NO_TITLE;
  else
    titleType_ = SINGLE;
  // Read one frame, see if the header is repeated.
  closeTraj();
  int nframes = TRAJIN_UNK;
  if (openTrajin()) return TRAJIN_ERR;
  if (titleType_ != NO_TITLE) infile_.Line();
  for (int at = 0; at != trajParm->Natom(); at++) {
    const char* atline = infile_.Line();
    if (atline == 0) {
      mprinterr("Error: Unexpected EOF when reading first frame of '%s'\n", atline);
      return TRAJIN_ERR;
    }
  }
  line2 = infile_.GetLine();
  if (!line2.empty()) {
    RemoveAllWhitespace( line2 );
    if (line2[0] == '#')
      titleType_ = MULTIPLE;
  } else
    nframes = 1;

  if (!line1.empty()) SetTitle( line1 );
  SetCoordInfo( CoordinateInfo(Box(), false, false, false) );
 
  set_ = 0;
 
  return nframes;
}

/** Read title. */
void Traj_XYZ::ReadTitle() {
  if (titleType_ == SINGLE) {
    infile_.Line();
    titleType_ = NO_TITLE;
  } else if (titleType_ == MULTIPLE)
    infile_.Line();
}

/** Read specified frame into given buffer. */
int Traj_XYZ::readXYZ(int set, int natom, double* xAddress) {
  // If an earlier set is being requested, reopen the file. 
  if (set < set_) {
    closeTraj();
    openTrajin();
  }
  // Perform any seeking needed
  while (set_ < set) {
    ReadTitle();
    for (int at = 0; at != natom; at++)
      infile_.Line();
    set_++;
  }
  ReadTitle(); 
  // Read coordinates into frame
  double* xyz = xAddress;
  for (int at = 0; at != natom; at++) {
    const char* ptr = infile_.Line();
    if (ptr == 0) return 1;
    if (sscanf(ptr, fmt_, xyz, xyz+1, xyz+2) != 3) return 1;
    xyz += 3;
  }
  set_++;
  return 0;
}

/** Read specified trajectory frame. */
int Traj_XYZ::readFrame(int set, Frame& frameIn) {
  return readXYZ(set, frameIn.Natom(), frameIn.xAddress());
}

/** Read velocities from specified frame. */
int Traj_XYZ::readVelocity(int set, Frame& frameIn) {
  return readXYZ(set, frameIn.Natom(), frameIn.vAddress());
}

/** Read forces from specified frame. */
int Traj_XYZ::readForce(int set, Frame& frameIn) {
  return readXYZ(set, frameIn.Natom(), frameIn.fAddress());
}

// -----------------------------------------------------------------------------
/** Write help. */
void Traj_XYZ::WriteHelp() {

}

/** Process write arguments. */
int Traj_XYZ::processWriteArgs(ArgList& argIn) {

  return 0;
}

/** Set up trajectory for write. */
int Traj_XYZ::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{
  titleType_ = SINGLE;
  ftype_ = ATOM_XYZ;
  return outfile_.OpenWrite( fname );
}

/** Write specified trajectory frame. */
int Traj_XYZ::writeFrame(int set, Frame const& frameOut) {
  if (titleType_ == SINGLE) {
    outfile_.Printf("#%s\n", Title().c_str());
    titleType_ = NO_TITLE;
  } else if (titleType_ == MULTIPLE)
    outfile_.Printf("#%s\n", Title().c_str());

  const double* xyz = frameOut.xAddress();
  if (ftype_ == ATOM_XYZ) {
    for (int at = 0; at != frameOut.Natom(); at++, xyz += 3)
      outfile_.Printf("%i %f %f %f\n", at+1, xyz[0], xyz[1], xyz[2]);
  } else if (ftype_ == XYZ) {
    for (int at = 0; at != frameOut.Natom(); at++, xyz += 3)
      outfile_.Printf("%f %f %f\n", xyz[0], xyz[1], xyz[2]);
  }
  return 0;
}

// =============================================================================
#ifdef MPI
/** Open trajectory for reading in parallel. */
int Traj_XYZ::parallelOpenTrajin(Parallel::Comm const& commIn) {
  return 1;
}

/** Open trajectory for writing in parallel. */
int Traj_XYZ::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return 1;
}

/** Set up trajectory for write in parallel. */
int Traj_XYZ::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{

  return 1;
}

/** Read frame in parallel. */
int Traj_XYZ::parallelReadFrame(int set, Frame& frameIn) {

  return 1;
}

/** Write frame in parallel. */
int Traj_XYZ::parallelWriteFrame(int set, Frame const& frameOut) {

  return 1;
}

/** Close trajectory in parallel. */
void Traj_XYZ::parallelCloseTraj() {

}
#endif
