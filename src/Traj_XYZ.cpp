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
  fmt_(0),
  hasBox_(false)
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

/** Read box information from given line. */
int Traj_XYZ::parseBoxLine(Box& xyzbox, ArgList const& boxline)
{
  double boxcrd[9];
  int iarg = 0;
  while (iarg < boxline.Nargs()) {
    int incr = 1;
    if (boxline[iarg] == "X:") {
      boxcrd[0] = convertToDouble(boxline[iarg+1]);
      boxcrd[1] = convertToDouble(boxline[iarg+2]);
      boxcrd[2] = convertToDouble(boxline[iarg+3]);
      incr = 4;
    } else if (boxline[iarg] == "Y:") {
      boxcrd[3] = convertToDouble(boxline[iarg+1]);
      boxcrd[4] = convertToDouble(boxline[iarg+2]);
      boxcrd[5] = convertToDouble(boxline[iarg+3]);
      incr = 4;
    } else if (boxline[iarg] == "Z:") {
      boxcrd[6] = convertToDouble(boxline[iarg+1]);
      boxcrd[7] = convertToDouble(boxline[iarg+2]);
      boxcrd[8] = convertToDouble(boxline[iarg+3]);
      incr = 4;
    }
    iarg += incr;
  }
  xyzbox.SetupFromUcell(boxcrd);
  return 0;
}

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
  Box xyzbox;
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
    // Determine if the comment line contains box information
    // 0    1  2   3  4      5     6
    // Conf 1. Box X: 19.660 0.000 0.000 Y: 0.000 19.660 0.000 Z: 0.000 0.000 19.660
    ArgList boxline(line2, " \r\n\t");
    if (boxline.Nargs() > 12) {
      if (boxline.hasKey("Box")) {
        mprintf("\tComment line appears to have box information.\n");
        hasBox_ = true;
        parseBoxLine(xyzbox, boxline);
      }
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
  SetCoordInfo( CoordinateInfo(xyzbox, false, false, false) );
 
  set_ = 0;
 
  return nframes;
}

/** Read title. */
void Traj_XYZ::ReadTitle(Box& xyzbox) {
  if (titleType_ == SINGLE) {
    file_.Line();
    titleType_ = NO_TITLE;
  } else if (titleType_ == MULTIPLE) {
    file_.Line();
  } else if (titleType_ == NATOM_COMMENT) {
    const char* ptr = file_.Line(); // #atoms
    if (ptr && hasBox_) {
      std::string line2 = file_.GetLine(); // comment
      ArgList boxline(line2, " \r\n\t");
      parseBoxLine(xyzbox, boxline);
    } else
      file_.Line(); // comment
  }
}

/** Read title. Skip box info */
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
int Traj_XYZ::readXYZ(int set, int natom, double* xAddress, Box& xyzbox) {
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
  ReadTitle(xyzbox); 
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
  return readXYZ(set, frameIn.Natom(), frameIn.xAddress(), frameIn.ModifyBox());
}

/** Read velocities from specified frame. */
int Traj_XYZ::readVelocity(int set, Frame& frameIn) {
  return readXYZ(set, frameIn.Natom(), frameIn.vAddress(), frameIn.ModifyBox());
}

/** Read forces from specified frame. */
int Traj_XYZ::readForce(int set, Frame& frameIn) {
  return readXYZ(set, frameIn.Natom(), frameIn.fAddress(), frameIn.ModifyBox());
}

// -----------------------------------------------------------------------------
/** Write help. */
void Traj_XYZ::WriteHelp() {
  mprintf("\tftype {namexyz|atomxyz|xyz}      : Choose either 'NAME X Y Z', 'ATOM X Y Z' (default) or 'X Y Z' output format.\n"
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
    else if (arg == "namexyz")
      ftype_ = NAME_XYZ;
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

  if (ftype_ == NAME_XYZ) {
    if (width_ == 0)
      width_ = 14;
    if (prec_ == -1)
      prec_ = 8;
  }
  TextFormat ffmt(TextFormat::DOUBLE, width_, prec_, 3);

  if (cInfoIn.HasBox() && ftype_ != NAME_XYZ) {
    mprintf("Warning: Box coordinates present but not using 'namexyz' format.\n"
            "Warning: Box coordinates will not be written.\n");
  }

  if (ftype_ == ATOM_XYZ)
    ofmt_ = "%i " + ffmt.Fmt();
  else if (ftype_ == XYZ)
    ofmt_ = ffmt.Fmt();
  else if (ftype_ == NAME_XYZ) {
    titleType_ = NATOM_COMMENT;
    ofmt_ = "%-7s " + ffmt.Fmt();
    names_.clear();
    names_.reserve( trajParm->Natom() );
    for (Topology::atom_iterator at = trajParm->begin(); at != trajParm->end(); ++at)
      names_.push_back( std::string(at->ElementName()) );
  }
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
  } else if (titleType_ == MULTIPLE) {
    file_.Printf("#%s\n", Title().c_str());
  } else if (titleType_ == NATOM_COMMENT) {
    file_.Printf("%i\n", frameOut.Natom());
    file_.Printf("Conf %i.", set+1);
    if (frameOut.BoxCrd().HasBox()) {
      Matrix_3x3 const& ucell = frameOut.BoxCrd().UnitCell();
      file_.Printf(" Box X: %.3f %.3f %.3f Y: %.3f %.3f %.3f Z: %.3f %.3f %.3f",
                   ucell[0], ucell[1], ucell[2],
                   ucell[3], ucell[4], ucell[5],
                   ucell[6], ucell[7], ucell[8]);
     }
    file_.Printf("\n");
  }

  const double* xyz = frameOut.xAddress();
  if (ftype_ == ATOM_XYZ) {
    for (int at = 0; at != frameOut.Natom(); at++, xyz += 3)
      file_.Printf(ofmt_.c_str(), at+1, xyz[0], xyz[1], xyz[2]);
  } else if (ftype_ == XYZ) {
    for (int at = 0; at != frameOut.Natom(); at++, xyz += 3)
      file_.Printf(ofmt_.c_str(), xyz[0], xyz[1], xyz[2]);
  } else if (ftype_ == NAME_XYZ) {
    for (int at = 0; at != frameOut.Natom(); at++, xyz += 3)
      file_.Printf(ofmt_.c_str(), names_[at].c_str(), xyz[0], xyz[1], xyz[2]);
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
