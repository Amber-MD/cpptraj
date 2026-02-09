#include "CMAP.h"
#include "../CpptrajStdio.h" // DEBUG
#include "../Frame.h"
#include "../ParameterTypes.h"
#include "../TorsionRoutines.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
CMAP::CMAP() {}

/// DEBUG - print gradients
static inline void print_grad(Vec3 const& dA, Vec3 const& dB, Vec3 const& dC, Vec3 const& dD)
{
  for (int i = 0; i < 3; i++)
    mprintf("DEBUG:\t\tdX %12.4f%12.4f%12.4f%12.4f\n", dA[i], dB[i], dC[i], dD[i]);
}

/** Calculate CMAP energy */
double CMAP::Ene_CMAP(CmapArray const& Cmaps, Frame const& frameIn)
const
{
  double ene_cmap = 0.0;

  for (CmapArray::const_iterator cmap = Cmaps.begin();
                                 cmap != Cmaps.end(); ++cmap)
  {
    const double* ixyz = frameIn.XYZ( cmap->A1() );
    const double* jxyz = frameIn.XYZ( cmap->A2() );
    const double* kxyz = frameIn.XYZ( cmap->A3() );
    const double* lxyz = frameIn.XYZ( cmap->A4() );
    const double* mxyz = frameIn.XYZ( cmap->A5() );

    // Calculate the dihedral angle (phi) and the derivatives of the
    // four coordinates with respect to phi. Remember this subroutine is
    // operating in radians.
    Vec3 dAphi, dBphi, dCphi, dDphi;
    double cosphi_ijkl, sinphi_ijkl;
    Torsion_and_part_deriv( ixyz, jxyz, kxyz, lxyz,
                            dAphi, dBphi, dCphi, dDphi,
                            cosphi_ijkl, sinphi_ijkl );
    print_grad(dAphi, dBphi, dCphi, dDphi);
    double phi = copysign(acos(cosphi_ijkl),sinphi_ijkl) * Constants::RADDEG;
    mprintf("DEBUG: Dihedral 1 %i %i %i %i = %g deg.\n", cmap->A1()+1, cmap->A2()+1, cmap->A3()+1, cmap->A4()+1, phi);

    // Calculate the dihedral angle (psi) and the derivatives of the
    // four coordinates with respect to psi. Remember this subroutine is
    // operating in radians.
    Vec3 dApsi, dBpsi, dCpsi, dDpsi;
    double cospsi_jklm, sinpsi_jklm;
    Torsion_and_part_deriv( jxyz, kxyz, lxyz, mxyz,
                            dApsi, dBpsi, dCpsi, dDpsi,
                            cospsi_jklm, sinpsi_jklm );
    print_grad(dApsi, dBpsi, dCpsi, dDpsi);
    double psi = copysign(acos(cospsi_jklm),sinpsi_jklm) * Constants::RADDEG;
    mprintf("DEBUG: Dihedral 2 %i %i %i %i = %g deg.\n", cmap->A2()+1, cmap->A3()+1, cmap->A4()+1, cmap->A5()+1, psi);
  }
  return ene_cmap;
}

/** This is reduced version of a typical cubic spline routine
  * that one may find in a recipe book of a numerical flavour.
  *
  * It is "reduced" since it is used here to obtain smooth derivatives
  * ONLY at discrete points on the CMAP grid, hence it never interpolates
  * between the points, therefore the coefficient a is always 1.0 and b is
  * always 0.0. In addition,  there is never a need to return the
  * interpolated value since it will be aways the same as the value passed to it.
  */
double CMAP::evaluate_cubic_spline(int step_size, std::vector<double> const& y,
                                   std::vector<double> const& y2,
                                   int grid_point)
{
  //Work out nearest complete grid point on the CMAP grid from xin
  //lo =  int( (xin - gridOrigin)/(step_size) ) + 1

  int lo = grid_point;
  //!write(6,'(a,I4)'),"Lo is: ",lo

  //!b = ( xin - ( (lo-1)*step_size + gridOrigin)  )/step_size
  //!a = 1-b

  //!write(6,'(a,f15.6)'),"a is : ",a
  //!write(6,'(a,f15.6)'),"b is : ",b

  double a = 1.0;
  double b = 0.0;
  //DEBUG

  //yout =   a*y(lo)                                        &
  //       + b*y(lo+1)                                      &
  //!      + (1/6)*(a*a*a-a)*(step_size*step_size)*y2(lo)   &
  //!      + (1/6)*(b*b*b-b)*(step_size*step_size)*y2(lo+1)

  double dyout =  (y[lo+1]-y[lo])/step_size
                - ((3*a*a-1)/6)*step_size*y2[lo]
                + ((3*b*b-1)/6)*step_size*y2[lo+1];

  //DEBUG
  //write(6,'(a,f15.6)'),"y(lo) is :", y(lo)
  //write(6,'(a,f15.6)'),"y(lo+1) is :", y(lo+1)
  //write(6,'(a,f15.6)'),"y2(lo) is :", y2(lo)
  //write(6,'(a,f15.6)'),"y2(lo+1) is :", y2(lo+1)
  //write(6,'(a,f15.6)'),"yout is :", yout
  //write(6,'(a)'),""
  return dyout;
}
