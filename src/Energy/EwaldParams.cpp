#include "EwaldParams.h"
#include "../Box.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Topology.h"

using namespace Cpptraj::Energy;

const double EwaldParams::INVSQRTPI_ = 1.0 / sqrt(Constants::PI);

/** CONSTRUCTOR */
EwaldParams::EwaldParams() :
  ew_coeff_(0.0),
  switch_width_(0.0),
  cutoff_(0.0),
  cut2_(0.0),
  cut2_0_(0.0),
  dsumTol_(0.0),
  debug_(0),
  NB_(0),
  sumq_(0),
  sumq2_(0)
{}


/** Determine Ewald coefficient from cutoff and direct sum tolerance.
  * Original Code: SANDER: findewaldcof
  */
double EwaldParams::FindEwaldCoefficient(double cutoff, double dsum_tol)
{
  // First get direct sum tolerance. How big must the Ewald coefficient be to
  // get terms outside the cutoff below tolerance?
  double xval = 0.5;
  int nloop = 0;
  double term = 0.0;
  do {
    xval = 2.0 * xval;
    nloop++;
    double yval = xval * cutoff;
    term = ErfcFxn::erfc_func(yval) / cutoff;
  } while (term >= dsum_tol);

  // Binary search tolerance is 2^-50
  int ntimes = nloop + 50;
  double xlo = 0.0;
  double xhi = xval;
  for (int i = 0; i != ntimes; i++) {
    xval = (xlo + xhi) / 2.0;
    double yval = xval * cutoff;
    double term = ErfcFxn::erfc_func(yval) / cutoff;
    if (term >= dsum_tol)
      xlo = xval;
    else
      xhi = xval;
  }
  mprintf("\t  Ewald coefficient determined using cut=%g, direct sum tol=%g is %g\n",
          cutoff, dsum_tol, xval);
  return xval;
}

/** Check some common input. */
int EwaldParams::CheckInput(Box const& boxIn, int debugIn, double cutoffIn, double dsumTolIn,
                      double ew_coeffIn, double switch_widthIn, double erfcTableDxIn, double skinnbIn)
{
  debug_ = debugIn;
  cutoff_ = cutoffIn;
  dsumTol_ = dsumTolIn;
  ew_coeff_ = ew_coeffIn;
  switch_width_ = switch_widthIn;
  double erfcTableDx = erfcTableDxIn;
  // Check input
  if (cutoff_ < Constants::SMALL) {
    mprinterr("Error: Direct space cutoff (%g) is too small.\n", cutoff_);
    return 1;
  }
  char dir[3] = {'X', 'Y', 'Z'};
  // NOTE: First 3 box parameters are X Y Z
  for (int i = 0; i < 3; i++) {
    if (cutoff_ > boxIn.Param((Box::ParamType)i)/2.0) {
      mprinterr("Error: Cutoff must be less than half the box length (%g > %g, %c)\n",
                cutoff_, boxIn.Param((Box::ParamType)i)/2.0, dir[i]);
      return 1;
    }
  }
  if (skinnbIn < 0.0) {
    mprinterr("Error: skinnb is less than 0.0\n");
    return 1;
  }
  if (switch_width_ < 0.0) switch_width_ = 0.0;
  if (switch_width_ > cutoff_) {
    mprinterr("Error: Switch width must be less than the cutoff.\n");
    return 1;
  }

  // Set defaults if necessary
  if (dsumTol_ < Constants::SMALL)
    dsumTol_ = 1E-5;
  if (DABS(ew_coeff_) < Constants::SMALL)
    ew_coeff_ = FindEwaldCoefficient( cutoff_, dsumTol_ );
  if (erfcTableDx <= 0.0) erfcTableDx = 1.0 / 5000;
  // TODO make this optional
  if (erfc_.FillErfcTable( erfcTableDx, 0.0, cutoff_*ew_coeff_*1.5 )) {
    mprinterr("Error: Could not set up spline table for ERFC\n");
    return 1;
  }
  
  // Calculate some common factors.
  cut2_ = cutoff_ * cutoff_;
  double cut0 = cutoff_ - switch_width_;
  cut2_0_ = cut0 * cut0;

  return 0;
}

/** Initialize */
int EwaldParams::InitEwald(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  if (CheckInput(boxIn, debugIn, pmeOpts.Cutoff(), pmeOpts.DsumTol(), pmeOpts.EwCoeff(),
                 pmeOpts.LJ_SwWidth(), pmeOpts.ErfcDx(), pmeOpts.SkinNB()))
    return 1;
  mprintf("\t  Cutoff= %g   Direct Sum Tol= %g   Ewald coeff.= %g  NB skin= %g\n",
          Cutoff(), DirectSumTol(), EwaldCoeff(), pmeOpts.SkinNB());
  if (LJ_SwitchWidth() > 0.0)
    mprintf("\t  LJ switch width= %g\n", LJ_SwitchWidth());
   //mprintf("\t  Erfc table dx= %g, size= %zu\n", erfcTableDx_, erfc_table_.size()/4);
  return 0;
}

/** Reserve space for selected atoms */
void EwaldParams::reserveRecipCoords(AtomMask const& maskIn) {
  coordsD_.reserve( maskIn.Nselected()*3 );
}

/** Fill recip coords with XYZ coords of selected atoms. */
void EwaldParams::FillRecipCoords(Frame const& frameIn, AtomMask const& maskIn)
{
  coordsD_.clear();
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm) {
    const double* XYZ = frameIn.XYZ( *atm );
    coordsD_.push_back( XYZ[0] );
    coordsD_.push_back( XYZ[1] );
    coordsD_.push_back( XYZ[2] );
  }
}

/** Convert charges to Amber units. Calculate sum of charges and squared charges.
  * Store LJ type indices for selected atoms.
  */
int EwaldParams::SetupEwald(Topology const& topIn, AtomMask const& maskIn) {
  NB_ = static_cast<NonbondParmType const*>( &(topIn.Nonbond()) );

  sumq_ = 0.0;
  sumq2_ = 0.0;
  Charge_.clear();
  Charge_.reserve( maskIn.Nselected() );
  TypeIndices_.clear();
  TypeIndices_.reserve( maskIn.Nselected() );
  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom) {
    double qi = topIn[*atom].Charge() * Constants::ELECTOAMBER;
    Charge_.push_back(qi);
    sumq_ += qi;
    sumq2_ += (qi * qi);
    // Store atom type indices for selected atoms.
    TypeIndices_.push_back( topIn[*atom].TypeIndex() );
  }
  //mprintf("DEBUG: sumq= %20.10f   sumq2= %20.10f\n", sumq_, sumq2_);
  //Setup_VDW_Correction( topIn, maskIn );
  reserveRecipCoords(maskIn);
  return 0;
}

/** Electrostatic self energy. This is the cancelling Gaussian plus the "neutralizing plasma". */
double EwaldParams::SelfEnergy(double volume) const {
//  t_self_.Start();
  double d0 = -ew_coeff_ * INVSQRTPI_;
  double ene = sumq2_ * d0;
//  mprintf("DEBUG: d0= %20.10f   ene= %20.10f\n", d0, ene);
  double factor = Constants::PI / (ew_coeff_ * ew_coeff_ * volume);
  double ee_plasma = -0.5 * factor * sumq_ * sumq_;
  ene += ee_plasma;
//  t_self_.Stop();
  return ene;
}

/** Electrostatic self energy on each atom. */
double EwaldParams::DecomposedSelfEnergy(Darray& atom_self, double volume) const {
  double d0 = -ew_coeff_ * INVSQRTPI_;
  double ene = sumq2_ * d0;
//  mprintf("DEBUG: d0= %20.10f   ene= %20.10f\n", d0, ene);
  double factor = Constants::PI / (ew_coeff_ * ew_coeff_ * volume);
  double ee_plasma = -0.5 * factor * sumq_ * sumq_;

  atom_self.resize( Charge_.size() );

  // Distribute the "neutrilizing plasma" to atoms equally
  for (unsigned int i = 0; i < Charge_.size(); i++)
  {
    atom_self[i] = Charge_[i]*Charge_[i]*d0 + ee_plasma/Charge_.size();
    //mprintf("DEBUG: Self energy atom %i = %f\n", i+1, atom_self[i]);
  }

  ene += ee_plasma;
  return ene;
}
