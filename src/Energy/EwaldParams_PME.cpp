#include "EwaldParams_PME.h"
#include "../AtomMask.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Frame.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
EwaldParams_PME::EwaldParams_PME() :
  order_(6)
{
  nfft_[0] = -1;
  nfft_[1] = -1;
  nfft_[2] = -1;
}

/** Set up PME parameters. */
int EwaldParams_PME::InitEwald(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  // Sanity check
  if (pmeOpts.Type() == EwaldOptions::REG_EWALD) {
    mprinterr("Internal Error: Options were set up for regular Ewald only.\n");
    return 1;
  }
  if (CheckInput(boxIn, debugIn, pmeOpts.Cutoff(), pmeOpts.DsumTol(), pmeOpts.EwCoeff(),
                 pmeOpts.LJ_SwWidth(), pmeOpts.ErfcDx(), pmeOpts.SkinNB()))
    return 1;
  nfft_[0] = pmeOpts.Nfft1();
  nfft_[1] = pmeOpts.Nfft2();
  nfft_[2] = pmeOpts.Nfft3();
  order_ = pmeOpts.SplineOrder();

  // Set defaults if necessary
  if (order_ < 1) order_ = 6;

  mprintf("\tParticle Mesh Ewald params:\n");
  mprintf("\t  Cutoff= %g   Direct Sum Tol= %g   Ewald coeff.= %g  NB skin= %g\n",
          Cutoff(), DirectSumTol(), EwaldCoeff(), pmeOpts.SkinNB());
    if (LJ_SwitchWidth() > 0.0)
    mprintf("\t  LJ switch width= %g\n", LJ_SwitchWidth());
  mprintf("\t  Bspline order= %i\n", order_);
  //mprintf("\t  Erfc table dx= %g, size= %zu\n", erfcTableDx_, erfc_table_.size()/4);
  mprintf("\t ");
  for (int i = 0; i != 3; i++)
    if (nfft_[i] == -1)
      mprintf(" NFFT%i=auto", i+1);
    else
      mprintf(" NFFT%i=%i", i+1, nfft_[i]);
  mprintf("\n");

  // Set up pair list
  //if (Setup_Pairlist(boxIn, pmeOpts.SkinNB())) return 1;

  return 0;
}

/** Setup PME calculation. */
int EwaldParams_PME::SetupEwald(Topology const& topIn, AtomMask const& maskIn) {
  CalculateCharges(topIn, maskIn);
  // Sanity check
  if (Charge().size() != (unsigned int)maskIn.Nselected()) {
    mprinterr("Internal Error: EwaldParams_PME::SetupEwald(): Charge/coords size mismatch (%zu vs %i\n",
              Charge().size(), maskIn.Nselected());
    return 1;
  }

  coordsD_.clear();
  coordsD_.reserve( maskIn.Nselected()*3 );
  return 0;
}

/** Put selected coords in separate array for recip calc. */
void EwaldParams_PME::FillRecipCoords(Frame const& frameIn, AtomMask const& maskIn)
{
  coordsD_.clear();
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm) {
    const double* XYZ = frameIn.XYZ( *atm );
    coordsD_.push_back( XYZ[0] );
    coordsD_.push_back( XYZ[1] );
    coordsD_.push_back( XYZ[2] );
  }
}