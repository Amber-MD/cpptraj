#include "CoordCovarMatrix.h"
#include "Atom.h"
#include "AtomMask.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
CoordCovarMatrix::CoordCovarMatrix() :
  nframes_(0),
  useMass_(false)
{}

/** Clear the matrix */
void CoordCovarMatrix::Clear() {
  covarMatrix_.clear();
  vect_.clear();
  mass_.clear();
  nframes_ = 0;
  clearMat();
}

/** Set up array sizes and masses. */
int CoordCovarMatrix::setupMat(std::vector<Atom> const& atoms,
                               AtomMask const& maskIn, bool useMassIn)
{
  useMass_ = useMassIn;
  // TODO more size error checking
  nframes_ = 0;
  vect_.assign(maskIn.Nselected(), Vec3(0.0));
  //Varray vect2(mask1_.Nselected(), Vec3(0.0));
  // Masses
  mass_.clear();
  if (useMassIn) {
    for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at)
      mass_.push_back( atoms[*at].Mass() );
  } else {
    for (int idx = 0; idx < maskIn.Nselected(); idx++)
      mass_.push_back( 1.0 );
  }
  
  return 0;
}

/** Debug print to file */
void CoordCovarMatrix::DebugPrint(const char* desc, CpptrajFile& outfile) const {
  if (desc != 0)
    outfile.Printf("DEBUG: CoordCovarMatrix: %s\n", desc);
  for (unsigned int row = 0; row < covarMatrix_.Nrows(); row++) {
    for (unsigned int col = 0; col < covarMatrix_.Ncols(); col++) {
      outfile.Printf(" %6.3f", covarMatrix_.element(col, row));
    }
    outfile.Printf("\n");
  }
}
