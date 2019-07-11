#include <cstdio> // sscanf
#include "Traj_CharmmCor.h"
#include "Topology.h"
#include "ArgList.h"
#include "Frame.h"
#include "CpptrajFile.h"
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
  mprintf("\tCOOR file: %i atoms\n", corAtom_);
  if (corAtom_ < 1) {
    mprinterr("Error: No atoms in CHARMM COOR file.\n");
    return TRAJIN_ERR;
  }
  if (corAtom_ > 99999)
    extendedFmt_ = true;
  if (corAtom_ != trajParm->Natom()) {
    mprinterr("Error: COOR file has %i atoms, associated topology '%s' has %i\n",
              corAtom_, trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  if (extendedFmt_)
    mprintf("\tCOOR file: extended format.\n");
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
      mprinterr("Error: Reading COOR atom %i\n", at+1);
      return 1;
    }
    int ncrd;
    if (extendedFmt_)
      ncrd = sscanf(buffer, "%*10i%*10i%*10s%*10s%20lf%20lf%20lf", xptr, xptr+1, xptr+2);
    else
      ncrd = sscanf(buffer, "%*5i%*5i%*5s%*5s%10lf%10lf%10lf", xptr, xptr+1, xptr+2);
    if (ncrd != 3) {
      mprinterr("Error: Reading coordinates for COOR atom %i (got %i)\n", at+1, ncrd);
      return 1;
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
// NOTE: The second-to-last field of the COOR format contains "residue ID",
//       which in CHARMM is written as a string but in practice is almost
//       always a number, so we will write it as such.
const char* Traj_CharmmCor::EXTENDED_FORMAT_ =
  "%10i%10i  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8i%20.10f\n";

const char* Traj_CharmmCor::REGULAR_FORMAT_ =
  "%5i%5i %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4i%10.5f\n";

void Traj_CharmmCor::WriteHelp() {
  mprintf("\tkeepext                : Keep filename extension; write '<name>.<num>.<ext>'\n"
          "\text                    : Use 'extended' format (default when > 99999 atoms.\n"
          "\tsegid <segid>          : Use <segid> as segment ID for all atoms.\n"
          "\tsegmask <mask> <segid> : Use <segid> as segment ID for atoms selected by <mask>.\n");
}

int Traj_CharmmCor::processWriteArgs(ArgList& argIn, DataSetList const& DSLin) {
  prependExt_ = argIn.hasKey("keepext");
  extendedFmt_ = argIn.hasKey("ext");
  std::string segid = argIn.GetStringKey("segid");
  if (!segid.empty()) {
    MaskSegPairs_.push_back("*");
    MaskSegPairs_.push_back( segid );
  }
  ArgList segpair = argIn.GetNstringKey("segmask", 2);
  while (segpair.Nargs() > 0) {
    if (segpair.Nargs() < 2) return 1;
    MaskSegPairs_.push_back( segpair[0] );
    MaskSegPairs_.push_back( segpair[1] );
    segpair = argIn.GetNstringKey("segmask", 2);
  }

  return 0;
}

int Traj_CharmmCor::setupTrajout(FileName const& fname, Topology* trajParm,
                               CoordinateInfo const& cInfoIn,
                               int NframesToWrite, bool append)
{
  if (trajParm==0) return 1;
  if (append) {
    mprinterr("Error: Append not supported for Charmm COOR.\n");
    return 1;
  }
  SetCoordInfo( cInfoIn );
  corTop_ = trajParm;
  corAtom_ = corTop_->Natom();
  if (corAtom_ > 99999 && !extendedFmt_) {
    mprintf("Info: > 99999 atoms; using extended COOR format.\n");
    extendedFmt_ = true;
  }
  if (extendedFmt_)
    outputFmt_ = EXTENDED_FORMAT_;
  else
    outputFmt_ = REGULAR_FORMAT_;
  // Set up file
  if (file_.SetupWrite( fname, debug_ )) return 1;
  if (NframesToWrite == 1)
    corWriteMode_ = SINGLE;
  else
    corWriteMode_ = MULTI;
  // Set up seg Ids.
  SegmentIds_.clear();
  SegmentIds_.reserve( trajParm->Nres() );
  if (!MaskSegPairs_.empty()) {
    // User-specified
    // First fill with default.
    SegmentIds_.assign(trajParm->Nres(), "PRO");
    for (Sarray::const_iterator sp = MaskSegPairs_.begin();
                                sp != MaskSegPairs_.end();
                                sp += 2)
    {
      // Mask, segid
      AtomMask mask( *sp );
      if (corTop_->SetupIntegerMask( mask )) return 1;
      mask.MaskInfo();
      for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at) {
        int res = (*corTop_)[ *at ].ResNum();
        SegmentIds_[res] = *(sp+1);
      }
    }
  } else if (trajParm->Res(0).ChainID() != ' ') {
    // Use chain IDs
    for (Topology::res_iterator res = trajParm->ResStart(); res != trajParm->ResEnd(); ++res)
      SegmentIds_.push_back( std::string(1, res->ChainID()) );
  } else
    // Default
    SegmentIds_.assign(trajParm->Nres(), "PRO");
  // Set default title if needed.
  if (Title().empty())
    SetTitle( "Cpptraj Generated COOR file.");
  else if (Title().size() > 78) {
    mprintf("Warning: Title for COOR too big, truncating.\n");
    std::string tmptitle = Title();
    tmptitle.resize(78);
    SetTitle( tmptitle );
  }

  return 0;
}
/*
void Traj_CharmmCor::writeTitle() {
  std::string buffer;
  const char* ptr = Title().c_str();
  unsigned int max = Title().size();
  unsigned int idx = 0;
  unsigned int col = 0;
  while (idx < max) {
    if (col == 0) {
      // Line start.
      file_.Printf("* ");
      col = 2;
    }
    // Find next newline or 80 chars, whatever comes first.
    unsigned int len = 0;
    for (; col < 80; col++, len++) {
      idx++;
      if ( idx == max || ptr[len] == '\n' ) {
        len++;
        break;
      }
*/

int Traj_CharmmCor::writeFrame(int set, Frame const& frameOut) {
  if (corWriteMode_ == MULTI) {
    if (file_.OpenWriteNumbered( set + 1, prependExt_ )) return 1;
  } else {
    if (file_.OpenFile()) return 1;
  }
  // Write title
  file_.Printf("* %s\n*\n", Title().c_str());
  // Write # atoms and optionally EXT
  if (extendedFmt_)
    file_.Printf("%10i  EXT\n", corTop_->Natom());
  else
    file_.Printf("%5i\n", corTop_->Natom());
  // Write each atom
  for (int res = 0; res < corTop_->Nres(); res++) {
    Residue const& Res = corTop_->Res(res);
    int resNum = Res.OriginalResNum();
    NameType const& resName = Res.Name();
    for (int at = Res.FirstAtom(); at != Res.LastAtom(); at++)
    {
      Atom const& Atm = (*corTop_)[at];
      const double* xyz = frameOut.XYZ( at );
      // TODO when is weight (WMAIN) non-zero?
      // Atom #, Res #,   Res name, At name, X, Y, Z,       segid,   res id, weight
      file_.Printf(outputFmt_, at+1, res+1, *resName, *(Atm.Name()),
                   xyz[0], xyz[1], xyz[2], SegmentIds_[res].c_str(), resNum, 0.0);
    }
  }
  file_.CloseFile();
  return 0;
}

// -----------------------------------------------------------------------------
void Traj_CharmmCor::Info() {
  mprintf("is a CHARMM COORdinates file");
  if (extendedFmt_) mprintf(" (extended format)");
}
