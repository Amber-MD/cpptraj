#include "EwaldParams_LJPME.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
EwaldParams_LJPME::EwaldParams_LJPME() :
   lw_coeff_(0.0)
{}

/** Set up LJPME parameters. */
int EwaldParams_LJPME::InitEwald(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  if (EwaldParams_PME::InitEwald(boxIn, pmeOpts, debugIn)) return 1;

  lw_coeff_ = pmeOpts.LwCoeff();

  // TODO do for C6 as well
  // TODO for C6 correction term
  if (lw_coeff_ < 0.0)
    lw_coeff_ = 0.0;
  else if (DABS(lw_coeff_) < Constants::SMALL)
    lw_coeff_ = EwaldCoeff();

  if (LJ_EwaldCoeff() > 0.0)
    mprintf("\t  LJ Ewald coeff.= %g\n", LJ_EwaldCoeff());
  return 0;
}
