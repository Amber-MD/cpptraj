#include <cmath> // sqrt
#include "Action_Projection.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Projection::Action_Projection() :
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

  // Check modes type
  if (modinfo_.Type() != DataSet_Matrix::COVAR &&
      modinfo_.Type() != DataSet_Matrix::MWCOVAR &&
      modinfo_.Type() != DataSet_Matrix::IDEA)
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
  for (int i = beg_; i < modinfo_.Nmodes() + beg_; ++i) {
    if (i == tgti) {
      --colwidth;
      if (colwidth < 5) colwidth = 5;
      tgti *= 10;
    }
    outfile_.Printf("%*s%i",colwidth, "Mode", i);
    if (modinfo_.Type() == DataSet_Matrix::IDEA)
      outfile_.Printf("                              ");
  }
  outfile_.Printf("\n");

  // Get mask
  mask_.SetMaskString( actionArgs.getNextMask() );

  mprintf("    PROJECTION: Calculating projection using modes %i to %i of file %s\n",
          beg_, end_, modinfo_.Legend().c_str());
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

// Action_Projection::setup()
int Action_Projection::setup() {
  // Setup mask
  if (currentParm->SetupIntegerMask( mask_ )) return 1;
  if (mask_.None()) {
    mprinterr("Error: projection: No atoms selected.\n");
    return 1;
  }
  mprintf("\t[%s] selected %i atoms.\n", mask_.MaskString(), mask_.Nselected());
  // Check # of selected atoms against modes info
  if ( modinfo_.Type() == DataSet_Matrix::COVAR || 
       modinfo_.Type() == DataSet_Matrix::MWCOVAR)
  {
    // Check if (3 * number of atoms in mask) and nvectelem agree
    int natom3 = mask_.Nselected() * 3;
    if ( natom3 != modinfo_.NavgCrd() ) {
      mprinterr("Error: projection: # selected coords != # avg coords in %s (%i)\n",
                natom3, modinfo_.Legend().c_str(), modinfo_.NavgCrd());
      return 1;
    }
    if ( natom3 != modinfo_.VectorSize() ) {
      mprinterr("Error: projection: # selected coords (%i) != eigenvector size (%i)\n",
                natom3, modinfo_.VectorSize() );
      return 1;
    }
  } else if ( modinfo_.Type() == DataSet_Matrix::IDEA ) {
    // Check if (number of atoms in mask) and nvectelem agree
    if (//mask_.Nselected() != modinfo_.Navgelem() ||
        mask_.Nselected() != modinfo_.VectorSize()) 
    {
      mprinterr("Error: projection: # selected atoms (%i) != eigenvector size (%i)\n",
                mask_.Nselected(), modinfo_.VectorSize() );
      return 1;
    }
  }

  // Precalc sqrt of mass for each coordinate
  sqrtmasses_.clear();
  if ( modinfo_.Type() == DataSet_Matrix::MWCOVAR ) {
    sqrtmasses_.reserve( mask_.Nselected() );
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      sqrtmasses_.push_back( sqrt( (*currentParm)[*atom].Mass() ) );
  } else {
    // If not MWCOVAR no mass-weighting necessary
    sqrtmasses_.resize( mask_.Nselected(), 1.0 );
  }

  return 0;
}

// Action_Projection::action()
int Action_Projection::action() {
    // If the current frame is less than start exit
  if (frameNum < start_) return 0;
  // If the current frame is greater than stop exit
  if (stop_!=-1 && frameNum >= stop_) return 0;
  // Update next target frame
  start_ += offset_;

  outfile_.Printf("%10i", frameNum+1);
  // Always start at first eigenvector element of first mode.
  const double* Vec = modinfo_.Eigenvector(0);
  // Project snapshots on modes
  if ( modinfo_.Type() == DataSet_Matrix::COVAR || 
       modinfo_.Type() == DataSet_Matrix::MWCOVAR ) 
  {
    for (int mode = 0; mode < modinfo_.Nmodes(); ++mode) {
      const double* Avg = modinfo_.AvgCrd();
      double proj = 0;
      std::vector<double>::const_iterator sqrtmass = sqrtmasses_.begin();
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      {
        const double* XYZ = currentFrame->XYZ( *atom );
        double mass = *(sqrtmass++);
        proj += (XYZ[0] - Avg[0]) * mass * Vec[0]; 
        proj += (XYZ[1] - Avg[1]) * mass * Vec[1]; 
        proj += (XYZ[2] - Avg[2]) * mass * Vec[2]; 
        Avg += 3;
        Vec += 3;
      }
      outfile_.Printf(" %9.3f", proj);
    }
  } else { // if modinfo_.Type() == IDEA
    for (int mode = 0; mode < modinfo_.Nmodes(); ++mode) {
      double proj1 = 0;
      double proj2 = 0;
      double proj3 = 0;
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      {
        const double* XYZ = currentFrame->XYZ(*atom);
        proj1 += XYZ[0] * *(Vec  );
        proj2 += XYZ[1] * *(Vec  );
        proj3 += XYZ[2] * *(Vec++);
      }
      outfile_.Printf(" %9.3f %9.3f %9.3f %9.3f", proj1, proj2, proj3,
                      sqrt(proj1*proj1 + proj2*proj2 + proj3*proj3) );
    }
  }
  outfile_.Printf("\n");
  return 0;
}
