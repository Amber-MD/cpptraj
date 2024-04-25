#include "CoordCovarMatrix.h"
#include "Atom.h"
#include "AtomMask.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "Frame.h"

/** CONSTRUCTOR - number of elements */
CoordCovarMatrix::CoordCovarMatrix() :
  nframes_(0),
  nelt_(0),
  useMass_(false)
{}

/** Clear the matrix */
void CoordCovarMatrix::Clear() {
  covarMatrix_.clear();
  nframes_ = 0;
  clearMat();
}

/** Setup mass array */
void CoordCovarMatrix::set_mass_array(Darray& mass, std::vector<Atom> const& atoms,
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

/** Debug print to file with default format. */
void CoordCovarMatrix::DebugPrint(const char* desc, CpptrajFile& outfile) const {
  DebugPrint(desc, outfile, " %6.3f");
}

/** Debug print to file */
void CoordCovarMatrix::DebugPrint(const char* desc, CpptrajFile& outfile, const char* fmt) const {
  if (desc != 0)
    outfile.Printf("DEBUG: CoordCovarMatrix: %s\n", desc);
  for (unsigned int row = 0; row < covarMatrix_.Nrows(); row++) {
    for (unsigned int col = 0; col < covarMatrix_.Ncols(); col++) {
      outfile.Printf(fmt, covarMatrix_.element(col, row));
    }
    outfile.Printf("\n");
  }
}

/** \return True if incoming array size is divisible by nelt_ */
bool CoordCovarMatrix::has_valid_size(Darray const& arrayIn) const {
  return ( (arrayIn.size() % nelt_) == 0);
}

/** Extract coordinates from given Frame into specified array. */
void CoordCovarMatrix::get_frame_coords(Darray& arrayIn, Frame const& frameIn, AtomMask const& maskIn)
const
{
  arrayIn.clear();
  unsigned int arraySize = (unsigned int)maskIn.Nselected() * nelt_;
  arrayIn.reserve( arraySize );
  for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at) {
    const double* XYZ = frameIn.XYZ(*at);
    arrayIn.push_back( XYZ[0] );
    arrayIn.push_back( XYZ[1] );
    arrayIn.push_back( XYZ[2] );
  }
}
