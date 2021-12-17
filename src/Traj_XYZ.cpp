#include <cstdio>
#include "Traj_XYZ.h"
#include "Topology.h"
#include "ArgList.h"
#include "Frame.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "TextFormat.h"

/// CONSTRUCTOR
Traj_XYZ::Traj_XYZ() :
  titleType_(UNKNOWN_TITLE),
  ftype_(UNKNOWN),
  set_(0),
  width_(0),
  prec_(-1),
  fmt_(0)
{}

/** Parse the given line, try to determine what format it corresponds to. */
Traj_XYZ::LineFmtType Traj_XYZ::DetermineLineFormat(std::string const& lineIn)
{
  if (lineIn.empty()) return UNKNOWN_LINE_FORMAT;
  ArgList line(lineIn, " \t\r\n");
  if (line.Nargs() == 1 && validInteger(line[0]))
    return SINGLE_INTEGER;
  else if (line.Nargs() == 3 && validDouble(line[0]) &&
                                validDouble(line[1]) &&
                                validDouble(line[2]))
    return THREE_DOUBLES;
  else if (line.Nargs() == 4 && validInteger(line[0]) &&
                                validDouble(line[1]) &&
                                validDouble(line[2]) &&
                                validDouble(line[3]))
    return INTEGER_AND_THREE_DOUBLES;
  else if (line.Nargs() == 4 && !validDouble(line[0]) && // TODO use isalpha(line[0][0])?
                                validDouble(line[1]) &&
                                validDouble(line[2]) &&
                                validDouble(line[3]))
    return STRING_AND_THREE_DOUBLES;
  else
    return UNKNOWN_LINE_FORMAT;
}

/** \param title Will be set to title line if title is present, cleared otherwise.
  * \param line1 First line.
  * \param line2 Second line.
  * \param line3 Third line.
  */
Traj_XYZ::Type Traj_XYZ::DetermineFormat(std::string& title,
                                         std::string const& line1,
                                         std::string const& line2,
                                         std::string const& line3)
{
  title.clear();
  LineFmtType l1fmt = DetermineLineFormat(line1);
  LineFmtType l2fmt = DetermineLineFormat(line2);
  LineFmtType l3fmt = DetermineLineFormat(line3);

  if (l1fmt == SINGLE_INTEGER && l2fmt == UNKNOWN_LINE_FORMAT && l3fmt == STRING_AND_THREE_DOUBLES)
  {
    title = line2;
    return NAME_XYZ;
  } else if (l1fmt == UNKNOWN_LINE_FORMAT) {
    if (l2fmt == INTEGER_AND_THREE_DOUBLES) {
      title = line1;
      return ATOM_XYZ;
    } else if (l2fmt == THREE_DOUBLES) {
      title = line1;
      return XYZ;
    }
  } else if (l1fmt == INTEGER_AND_THREE_DOUBLES) {
    return ATOM_XYZ;
  } else if (l1fmt == THREE_DOUBLES)
    return XYZ;
  
/*
  std::string line = line1;
  // This line can be a title line, atom XYZ, XYZ, or # atoms
  RemoveLeadingWhitespace( line );
  RemoveTrailingWhitespace( line );
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
  //mprintf("DEBUG: '%s' '%s' '%s' '%s'\n", tk0, tk1, tk2, tk3);
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
  }*/
  return UNKNOWN;
}

/** Identify trajectory format. File should be setup for READ */
bool Traj_XYZ::ID_TrajFormat(CpptrajFile& fileIn) {
  // File must already be set up for read.
  if (fileIn.OpenFile()) return false;
  std::string line1 = fileIn.GetLine();
  std::string line2 = fileIn.GetLine();
  std::string line3 = fileIn.GetLine();
  std::string title;
  if ( DetermineFormat(title, line1, line2, line3) == UNKNOWN ) return false;

  return true;
}

/** Print trajectory info to stdout. */
void Traj_XYZ::Info() {
  switch (ftype_) {
    case UNKNOWN  :
    case XYZ      : mprintf("is an XYZ trajectory"); break;
    case ATOM_XYZ : mprintf("is an Atom-XYZ trajectory"); break;
    case NAME_XYZ : mprintf("is a regular XYZ trajectory"); break;
  }
}

/** Close file. */
void Traj_XYZ::closeTraj() {
  file_.CloseFile();
}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int Traj_XYZ::openTrajin() {
  if (file_.OpenFile()) return 1;
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

const char* Traj_XYZ::FMT_NAME_XYZ_ = "%*s %lf %lf %lf";

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int Traj_XYZ::setupTrajin(FileName const& fname, Topology* trajParm)
{
  if (file_.OpenFileRead( fname )) return TRAJIN_ERR;
  // Initial format determiniation.
  std::string title;
  std::string line1 = file_.GetLine();
  std::string line2 = file_.GetLine();
  std::string line3 = file_.GetLine();
  ftype_ = DetermineFormat(title, line1, line2, line3);
  switch (ftype_) {
    case XYZ      :
      fmt_ = FMT_XYZ_;
      mprintf("\tCoordinate lines are <X> <Y> <Z>\n");
      break;
    case ATOM_XYZ :
      fmt_ = FMT_ATOM_XYZ_;
      mprintf("\tCoordinate lines are <#> <X> <Y> <Z>\n");
      break;
    case NAME_XYZ :
      fmt_ = FMT_NAME_XYZ_;
      mprintf("\tCoordinate lines are <name> <X> <Y> <Z>\n");
      break;
    case UNKNOWN  :
      mprinterr("Internal Error: '%s' does not appear to be XYZ format anymore.\n",
                file_.Filename().full());
      return TRAJIN_ERR;
  }
  // Initial header type determination
  if (title.empty())
    titleType_ = NO_TITLE;
  else if (ftype_ == NAME_XYZ) {
    titleType_ = NATOM_COMMENT;
    int nat = convertToInteger(line1);
    if (nat != trajParm->Natom()) {
      mprinterr("Error: # atoms in XYZ (%i) does not match # atoms in topology '%s' (%i)\n",
                nat, trajParm->c_str(), trajParm->Natom());
      return TRAJIN_ERR;
    }
  } else
    titleType_ = SINGLE;
  closeTraj();
  int nframes = TRAJIN_UNK;

  // Read one frame, see if the header is repeated.
  if (titleType_ == SINGLE) {
    if (openTrajin()) return TRAJIN_ERR;
    if (titleType_ == SINGLE)
      file_.Line();
    for (int at = 0; at != trajParm->Natom(); at++) {
      const char* atline = file_.Line();
      if (atline == 0) {
        mprinterr("Error: Unexpected EOF when reading first frame of '%s'\n", atline);
        return TRAJIN_ERR;
      }
    }
    line2 = file_.GetLine();
    if (!line2.empty()) {
      if (DetermineLineFormat(line2) == UNKNOWN_LINE_FORMAT)
        titleType_ = MULTIPLE;
      //RemoveLeadingWhitespace( line2 );
      //RemoveTrailingWhitespace( line2 );
      //if (line2[0] == '#')
      //  titleType_ = MULTIPLE;
    } else
      nframes = 1;
  }

  switch (titleType_) {
    case UNKNOWN_TITLE :
    case NO_TITLE      : mprintf("\tNo title detected.\n"); break;
    case SINGLE        : mprintf("\tSingle title detected.\n"); break;
    case MULTIPLE      : mprintf("\tTitle before each frame detected.\n"); break;
    case NATOM_COMMENT : mprintf("\tStandard # atoms followed by comment detected.\n"); break;
  }

  if (!title.empty()) SetTitle( title );
  SetCoordInfo( CoordinateInfo(Box(), false, false, false) );
 
  set_ = 0;
 
  return nframes;
}

/** Read title. */
void Traj_XYZ::ReadTitle() {
  if (titleType_ == SINGLE) {
    file_.Line();
    titleType_ = NO_TITLE;
  } else if (titleType_ == MULTIPLE) {
    file_.Line();
  } else if (titleType_ == NATOM_COMMENT) {
    file_.Line(); // #atoms
    file_.Line(); // comment
  }
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
      file_.Line();
    set_++;
  }
  ReadTitle(); 
  // Read coordinates into frame
  double* xyz = xAddress;
  for (int at = 0; at != natom; at++) {
    const char* ptr = file_.Line();
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
  mprintf("\tftype {atomxyz|xyz}              : Choose either 'ATOM X Y Z' (default) or 'X Y Z' output format.\n"
          "\ttitletype {none|single|perframe} : No title, one title (default), or title before every frame.\n"
          "\twidth <#>                        : Output format width.\n"
          "\tprec <#>                         : Output format precision.\n");
}

/** Process write arguments. */
int Traj_XYZ::processWriteArgs(ArgList& argIn, DataSetList const& DSLin) {
  std::string arg = argIn.GetStringKey("ftype");
  if (!arg.empty()) {
    if (arg == "xyz")
      ftype_ = XYZ;
    else if (arg == "atomxyz")
      ftype_ = ATOM_XYZ;
    else {
      mprinterr("Error: Unrecognized key for 'ftype' = '%s'\n", arg.c_str());
      return 1;
    }
  }
  arg = argIn.GetStringKey("titletype");
  if (!arg.empty()) {
    if (arg == "none")
      titleType_ = NO_TITLE;
    else if (arg == "single")
      titleType_ = SINGLE;
    else if (arg == "perframe")
      titleType_ = MULTIPLE;
    else {
      mprinterr("Error: Unrecognized key for 'titletype' = '%s'\b", arg.c_str());
      return 1;
    }
  }
  width_ = argIn.getKeyInt("width", width_);
  prec_ = argIn.getKeyInt("prec", prec_);
  return 0;
}

/** Set up trajectory for write. */
int Traj_XYZ::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{
  if (titleType_ == UNKNOWN_TITLE)
    titleType_ = SINGLE;
  if (ftype_ == UNKNOWN)
    ftype_ = ATOM_XYZ;
  TextFormat ffmt(TextFormat::DOUBLE, width_, prec_, 3);
  if (ftype_ == ATOM_XYZ)
    ofmt_ = "%i " + ffmt.Fmt();
  else if (ftype_ == XYZ)
    ofmt_ = ffmt.Fmt();
  ofmt_.append("\n");
  //mprintf("DEBUG: output format string: '%s'\n", ofmt_.c_str()); 
  if (titleType_ != NO_TITLE && Title().empty())
    SetTitle("Cpptraj Generated XYZ file.");

  return file_.OpenWrite( fname );
}

/** Write specified trajectory frame. */
int Traj_XYZ::writeFrame(int set, Frame const& frameOut) {
  if (titleType_ == SINGLE) {
    file_.Printf("#%s\n", Title().c_str());
    titleType_ = NO_TITLE;
  } else if (titleType_ == MULTIPLE)
    file_.Printf("#%s\n", Title().c_str());

  const double* xyz = frameOut.xAddress();
  if (ftype_ == ATOM_XYZ) {
    for (int at = 0; at != frameOut.Natom(); at++, xyz += 3)
      file_.Printf(ofmt_.c_str(), at+1, xyz[0], xyz[1], xyz[2]);
  } else if (ftype_ == XYZ) {
    for (int at = 0; at != frameOut.Natom(); at++, xyz += 3)
      file_.Printf(ofmt_.c_str(), xyz[0], xyz[1], xyz[2]);
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
