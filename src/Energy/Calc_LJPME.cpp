#include "Calc_LJPME.h"

using namespace Cpptraj::Energy;

Calc_LJPME::Calc_LJPME() :
  Recip_(PME_Recip::COULOMB),
  LJrecip_(PME_Recip::LJ)
{}

