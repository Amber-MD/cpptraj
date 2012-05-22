#include <cmath> // sqrt
#include "Action_Projection.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Projection::Action_Projection() :
  type_(ModesInfo::MT_UNKNOWN),
  beg_(0),
  end_(0),
  start_(0),
  stop_(-1),
  offset_(1)
{}

/** projection modes <modesfile> out <outfile>
  *            [beg <beg>] [end <end>] [<mask>]
  *            [start <start>] [stop <stop>] [offset <offset>]
  */
int Action_Projection::init() {
  // Get ibeg, iend, start, stop, offset
  beg_ = actionArgs.getKeyInt("beg", 1);
  end_ = actionArgs.getKeyInt("end", 2);
  start_ = actionArgs.getKeyInt("start", 1);
  // Frames start from 0
  --start_;
  stop_ = actionArgs.getKeyInt("stop", -1);
  offset_ = actionArgs.getKeyInt("offset", 1);

  // Get modes from file
  std::string modesname = actionArgs.GetStringKey("modes");
  if (modesname.empty()) {
    mprinterr("Error: projection: no modes file specified ('modes <filename>')\n");
    return 1;
  }
  if (modinfo_.ReadEvecFile( modesname, beg_, end_ )) return 1;
  // NOTE: This warning may be redundant
  if (modinfo_.Nvect() != (end_ - beg_ + 1)) {
    mprintf("Warning: projection: # read evecs is %i, # requested is %i\n",
            modinfo_.Nvect(), (end_ - beg_ + 1));
  }

  // Check modes type
  type_ = modinfo_.Mtype();
  if (type_ != ModesInfo::MT_COVAR &&
      type_ != ModesInfo::MT_MWCOVAR &&
      type_ != ModesInfo::MT_IDEA)
  {
    mprinterr("Error: projection: evecs type is not COVAR, MWCOVAR, or IDEA.\n");
    return 1;
  }

  // Output Filename
  std::string filename_ = actionArgs.GetStringKey("out");
  if (filename_.empty()) {
    mprinterr("Error: projection: No output file specified ('out <filename>')\n");
    return 1;
  }

  // Open output file, write header
  if (outfile_.OpenWrite( filename_ )) return 1;
  outfile_.Printf("Projection of snapshots onto modes\n%10s", "Snapshot");
  int colwidth = 9;
  int tgti = 10;
  for (int i = beg_; i < modinfo_.Nvect() + beg_; ++i) {
    if (i == tgti) {
      --colwidth;
      if (colwidth < 5) colwidth = 5;
      tgti *= 10;
    }
    outfile_.Printf("%*s%i",colwidth, "Mode", i);
    if (type_ == ModesInfo::MT_IDEA)
      outfile_.Printf("                              ");
  }
  outfile_.Printf("\n");

  // Get mask
  mask_.SetMaskString( actionArgs.getNextMask() );

  mprintf("    PROJECTION: Calculating projection using modes %i to %i of file %s\n",
          beg_, end_, modinfo_.c_str());
  mprintf("                Results are written to %s\n", filename_.c_str());
  mprintf("                Start: %i", start_+1);
  if (stop_!=-1)
    mprintf("  Stop: %i", stop_);
  else
    mprintf("  Stop: Last Frame");
  if (offset_!=1)
    mprintf("  Offset: %i", offset_);
  mprintf("\n");
  mprintf("                Atom Mask: [%s]\n", mask_.MaskString());

  return 0;
}

int Action_Projection::setup() {
  // Setup mask
  if (currentParm->SetupIntegerMask( mask_ )) return 1;
  if (mask_.None()) {
    mprinterr("Error: projection: No atoms selected.\n");
    return 1;
  }
  mprintf("\t[%s] selected %i atoms.\n", mask_.MaskString(), mask_.Nselected());
  // Check # of selected atoms against modes info
  if ( type_ == ModesInfo::MT_COVAR || type_ == ModesInfo::MT_MWCOVAR)
  {
    // Check if (3 * number of atoms in mask) and nvectelem agree
    int natom3 = mask_.Nselected() * 3;
    if ( natom3 != modinfo_.Navgelem() || natom3 != modinfo_.NvectElem() )
    {
      mprinterr("Error: projection: # selected coords (%i) != # evec elements (%i, %i)\n",
                natom3, modinfo_.Navgelem(), modinfo_.NvectElem() );
      return 1;
    }
  } else if ( type_ == ModesInfo::MT_IDEA ) {
    // Check if (number of atoms in mask) and nvectelem agree
    if (mask_.Nselected() != modinfo_.Navgelem() ||
        mask_.Nselected() != modinfo_.NvectElem()) 
    {
      mprinterr("Error: projection: # selected atoms (%i) != # evec elements (%i, %i)\n",
                mask_.Nselected(), modinfo_.Navgelem(), modinfo_.NvectElem() );
      return 1;
    }
  }

  // Precalc sqrt of mass for each coordinate
  sqrtmasses_.clear();
  if ( type_ == ModesInfo::MT_MWCOVAR ) {
    sqrtmasses_.reserve( mask_.Nselected() );
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      sqrtmasses_.push_back( sqrt( (*currentParm)[*atom].Mass() ) );
  } else {
    // If not MWCOVAR no mass-weighting necessary
    sqrtmasses_.resize( mask_.Nselected(), 1.0 );
  }

  return 0;
}

int Action_Projection::action() {
    // If the current frame is less than start exit
  if (frameNum < start_) return 0;
  // If the current frame is greater than stop exit
  if (stop_!=-1 && frameNum >= stop_) return 0;
  // Update next target frame
  start_ += offset_;

  outfile_.Printf("%10i", frameNum+1);
  // Project snapshots on modes
  if ( type_ == ModesInfo::MT_COVAR || type_ == ModesInfo::MT_MWCOVAR )
    modinfo_.ProjectCovar(outfile_, *currentFrame, mask_, sqrtmasses_);
  else // if type_ == MT_IDEA
    modinfo_.ProjectIDEA(outfile_, *currentFrame, mask_);
     
  return 0;
}
