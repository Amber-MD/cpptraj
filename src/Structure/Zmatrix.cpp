#include <vector>
#include "Zmatrix.h"
#include "InternalCoords.h"
#include "../Frame.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Zmatrix::Zmatrix() :
  seed0_(InternalCoords::NO_ATOM),
  seed1_(InternalCoords::NO_ATOM),
  seed2_(InternalCoords::NO_ATOM)
{}

/** Setup Zmatrix from coordinates/topology. */
int Zmatrix::SetFromFrame(Frame const& frameIn, Topology const& topIn)
{
  IC_.clear();
  // First seed is first atom. Store XYZ.
  IC_.push_back( InternalCoords( frameIn.XYZ(0) ) );
  seed0_ = 0;

  return 0;
}
