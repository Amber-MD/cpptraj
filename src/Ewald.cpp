#include "Ewald.h"
#include "CpptrajStdio.h"
#include <cmath>

Ewald::Ewald() : sumq_(0.0), sumq2_(0.0), ew_coeff_(0.0)
{}

// Original code: SANDER: erfcfun.F90
double Ewald::erfc_func(double xIn) const {
  double erfc;
  double absx;
  if (xIn < 0.0)
    absx = -xIn;
  else
    absx = xIn;
    
  if (xIn > 26.0)
    erfc = 0.0;
  else if (xIn < -5.5)
    erfc = 2.0;
  else if (absx <= 0.5) {
    double cval = xIn * xIn;
    double pval = ((-0.356098437018154E-1*cval+0.699638348861914E1)*cval + 0.219792616182942E2) *
                  cval + 0.242667955230532E3;
    double qval = ((cval+0.150827976304078E2)*cval+0.911649054045149E2)*cval + 0.215058875869861E3;
    double erf = xIn * pval/qval;
    erfc = 1.0 - erf;
  } else if (absx < 4.0) {
    double cval = absx;
    double pval=((((((-0.136864857382717E-6*cval+0.564195517478974)*cval+
                     0.721175825088309E1)*cval+0.431622272220567E2)*cval+
                   0.152989285046940E3)*cval+0.339320816734344E3)*cval+
                 0.451918953711873E3)*cval+0.300459261020162E3;
    double qval=((((((cval+0.127827273196294E2)*cval+0.770001529352295E2)*cval+
                    0.277585444743988E3)*cval+0.638980264465631E3)*cval+
                  0.931354094850610E3)*cval+0.790950925327898E3)*cval+
                0.300459260956983E3;
    double nonexperfc;
    if ( xIn > 0.0 )
      nonexperfc = pval/qval;
    else
      nonexperfc = 2.0*exp(xIn*xIn) - pval/qval;
    erfc = exp(-absx*absx)*nonexperfc;
  } else {
    double cval = 1.0/(xIn*xIn);
    double pval = (((0.223192459734185E-1*cval+0.278661308609648)*cval+
                    0.226956593539687)*cval+0.494730910623251E-1)*cval+
                  0.299610707703542E-2;
    double qval = (((cval+0.198733201817135E1)*cval+0.105167510706793E1)*cval+
                   0.191308926107830)*cval+0.106209230528468E-1;
    cval = (-cval*pval/qval + 0.564189583547756)/absx;
    double nonexperfc;
    if ( xIn > 0.0 )
      nonexperfc = cval;
    else
      nonexperfc = 2.0*exp(xIn*xIn) - cval;
    erfc = exp(-absx*absx)*nonexperfc;
  }
  return erfc;
}
  
// Original Code: SANDER: findewaldcof
void Ewald::FindEwaldCoefficient(double cutoff, double dsum_tol)
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
    term = erfc_func(yval) / cutoff;
  } while (term >= dsum_tol);

  // Binary search tolerance is 2^-50
  int ntimes = nloop + 50;
  double xlo = 0.0;
  double xhi = xval;
  for (int i = 0; i != ntimes; i++) {
    xval = (xlo + xhi) / 2.0;
    double yval = xval * cutoff;
    double term = erfc_func(yval) / cutoff;
    if (term >= dsum_tol)
      xlo = xval;
    else
      xhi = xval;
  }
  ew_coeff_ = xval;
  mprintf("DEBUG: Ewald coefficient for cut=%g, direct sum tol=%g is %g\n",
          cutoff, dsum_tol, ew_coeff_);
}
