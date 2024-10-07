#include "RecipCoords.h"
#include "../AtomMask.h"
#include "../Frame.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
RecipCoords::RecipCoords() {}

/** Reserve space for selected atoms */
int RecipCoords::ReserveForSelected(AtomMask const& maskIn) {
  coordsD_.reserve( maskIn.Nselected()*3 );
  return 0;
}

/** Fill recip coords with XYZ coords of selected atoms. */
void RecipCoords::FillRecipCoords(Frame const& frameIn, AtomMask const& maskIn)
{
  coordsD_.clear();
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm) {
    const double* XYZ = frameIn.XYZ( *atm );
    coordsD_.push_back( XYZ[0] );
    coordsD_.push_back( XYZ[1] );
    coordsD_.push_back( XYZ[2] );
  }
}
