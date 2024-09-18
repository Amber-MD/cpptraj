#include "EwaldParams.h"
#include "../Box.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Energy;

static inline double DABS(double xIn) { if (xIn < 0.0) return -xIn; else return xIn; }

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
  mprintf("\tEwald coefficient for cut=%g, direct sum tol=%g is %g\n",
          cutoff, dsum_tol, xval);
  return xval;
}

/** Check some common input. */
int EwaldParams::CheckInput(Box const& boxIn, int debugIn, double cutoffIn, double dsumTolIn,
                      double ew_coeffIn, double lw_coeffIn, double switch_widthIn,
                      double erfcTableDxIn, double skinnbIn)
{
  debug_ = debugIn;
  cutoff_ = cutoffIn;
  dsumTol_ = dsumTolIn;
  ew_coeff_ = ew_coeffIn;
  lw_coeff_ = lw_coeffIn;
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
  // TODO do for C6 as well
  // TODO for C6 correction term
  if (lw_coeff_ < 0.0)
    lw_coeff_ = 0.0;
  else if (DABS(lw_coeff_) < Constants::SMALL)
    lw_coeff_ = ew_coeff_;

  // Calculate some common factors.
  cut2_ = cutoff_ * cutoff_;
  double cut0 = cutoff_ - switch_width_;
  cut2_0_ = cut0 * cut0;

  return 0;
}

