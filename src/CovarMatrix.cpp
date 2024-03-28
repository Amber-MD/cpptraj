#include "CovarMatrix.h"
#include "Atom.h"
#include "AtomMask.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "DataSet_2D.h"

/** CONSTRUCTOR */
CovarMatrix::CovarMatrix(unsigned int n) :
  covarMatrix_(0),
  nelt_(n),
  useMass_(false)
{}

/** Clear the matrix */
void CovarMatrix::Clear() {
  if (covarMatrix_ != 0) covarMatrix_->Clear();
  nframes_ = 0;
  clearMat();
}

/** Setup mass array */
void CovarMatrix::set_mass_array(Darray& mass, std::vector<Atom> const& atoms,
                                 AtomMask const& maskIn, bool useMassIn)
{
  nframes_ = 0;
  useMass_ = useMassIn;
  mass.clear();
  mass.reserve( maskIn.Nselected() );
  if (useMassIn) {
    for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at)
      mass.push_back( atoms[*at].Mass() );
  } else {
    for (int idx = 0; idx < maskIn.Nselected(); idx++)
      mass.push_back( 1.0 );
  }
}

/** Debug print to file */
void CovarMatrix::DebugPrint(const char* desc, CpptrajFile& outfile) const {
  if (desc != 0)
    outfile.Printf("DEBUG: CoordCovarMatrix: %s\n", desc);
  if (covarMatrix_ == 0) {
    mprintf("Warning: No covar matrix allocated.\n");
    outfile.Printf("No covar matrix allocated.\n");
    return;
  }

  for (unsigned int row = 0; row < covarMatrix_->Nrows(); row++) {
    for (unsigned int col = 0; col < covarMatrix_->Ncols(); col++) {
      outfile.Printf(" %6.3f", covarMatrix_->GetElement(col, row));
    }
    outfile.Printf("\n");
  }
}
