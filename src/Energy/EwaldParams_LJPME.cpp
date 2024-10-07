#include "EwaldParams_LJPME.h"
#include "../AtomMask.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Topology.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
EwaldParams_LJPME::EwaldParams_LJPME() :
   lw_coeff_(0.0),
   ljpme_self_(0.0)
{}

/** Set up LJPME parameters. */
int EwaldParams_LJPME::InitEwald(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  if (EwaldParams::InitEwald(boxIn, pmeOpts, debugIn)) return 1;

  lw_coeff_ = pmeOpts.LwCoeff();

  // TODO do for C6 as well
  // TODO for C6 correction term
  if (lw_coeff_ < 0.0) {
    mprinterr("Internal Error: LJ PME requested but LJ Ewald coefficient is < 0\n");
    return 1;
  }
  //  lw_coeff_ = 0.0;
  //else
  if (DABS(lw_coeff_) < Constants::SMALL)
    lw_coeff_ = EwaldCoeff();

  if (LJ_EwaldCoeff() > 0.0)
    mprintf("\t  LJ Ewald coeff.= %g\n", LJ_EwaldCoeff());
  return 0;
}

/** Setup LJPME calculation. */
int EwaldParams_LJPME::SetupEwald(Topology const& topIn, AtomMask const& maskIn) {
  if (EwaldParams::SetupEwald(topIn, maskIn)) return 1;

  // Calcuate C6 parameters
  if (lw_coeff_ > 0.0) {
    for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
    {
      double rmin = topIn.GetVDWradius( *atom );
      double eps  = topIn.GetVDWdepth( *atom );
      Cparam_.push_back( 8.0 * (rmin*rmin*rmin) * sqrt(2 * eps) );
      if (Debug() > 0)
        mprintf("DEBUG: C6 param atom %8i = %16.8f\n", *atom+1, Cparam_.back());
    }
  } else {
    mprinterr("Internal Error: LJ PME setup called with LJ Ewald coefficient <= 0\n");
    return 1;
  }
  // Lennard-Jones self energy. (used to be Ewald::Self6())
  double ew2 = lw_coeff_ * lw_coeff_;
  double ew6 = ew2 * ew2 * ew2;
  double c6sum = 0.0;
  for (Darray::const_iterator it = Cparam_.begin(); it != Cparam_.end(); ++it)
    c6sum += ew6 * (*it * *it);
  ljpme_self_ = c6sum / 12.0;

  return 0;
}

/** Calculate decomposed self6 energies. */
void EwaldParams_LJPME::CalcDecomposedSelf6Energy() {
  atom_vdwself6_.clear();
  atom_vdwself6_.reserve( Cparam_.size() );
  double ew2 = lw_coeff_ * lw_coeff_;
  double ew6 = ew2 * ew2 * ew2;
  double c6sum = 0.0;
  for (Darray::const_iterator it = Cparam_.begin(); it != Cparam_.end(); ++it)
  {
    c6sum += ew6 * (*it * *it);
    atom_vdwself6_.push_back(ew6 * (*it * *it)/12.0);
  }
}
