#include <cmath>
#include "AtomType.h"
#include "Constants.h"

bool AtomType::operator==(AtomType const& rhs) const {
  return ( (fabs(radius_ - rhs.radius_) < Constants::SMALL) &&
           (fabs(depth_  - rhs.depth_ ) < Constants::SMALL) );
}

NonbondType AtomType::Combine_LB(AtomType const& rhs) const {
  double dR = radius_ + rhs.radius_;
  double dE = sqrt( depth_ * rhs.depth_ );
  double dR2 = dR * dR;
  double dR6 = dR2 * dR2 * dR2;
  double dER6 = dE * dR6;
  return NonbondType( dER6*dR6, 2.0*dER6 );
}
