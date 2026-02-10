#include "CMAP.h"
#include "../CpptrajStdio.h" // DEBUG
#include "../Frame.h"
#include "../ParameterTypes.h"
#include "../Topology.h"
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

/** Generates the set cubic splining coefficients from a 1D array with
  * the assumption that the distance between the points is constant. This
  * interpolation method ensures that the derivative of this spline is
  * continuous across the boundary of two intervals.
  *
  * These coefficients are only calculated *once*, but are later used on
  * multiple occassions by another subroutine to interpolate anywhere
  * between two points in this 1D array.
  */
void CMAP::generate_cubic_spline(int n, int step_size, std::vector<double> const& y, std::vector<double>& y2)
{
  // tmp array used internally for the decomposition loop
  std::vector<double> tmp(n);
  // Lower and upper boundaries are natural
  y2[0]  = 0.0;
  tmp[0] = 0.0;

  for (int i = 1; i < n-1; i++)
  {
    // y2 is used initially as a temp storage array
    double p = 0.5*y2[i-1] + 2.0;
    //Debug
    y2[i]  = -0.5/p;
    tmp[i] = ( y[i+1]-2.0*y[i] + y[i-1] ) / (double)step_size;
    tmp[i] = ( 3.0* (tmp[i]/(double)step_size) - 0.5 *tmp[ i-1 ] )/p;
    mprintf("%6i y2= %16.8f  tmp= %16.8f\n", i, y2[i], tmp[i]);
  }

  // Set the upper boundary
  y2[n-1] = 0.0;

  for (int i = n-2; i > -1; i--)
  {
    y2[i]=y2[i]*y2[i+1]+tmp[i];
    mprintf("%6i %16.8f\n", i, y2[i]);
    //mprintf("\t\t\tn is %i:\n",n);
    //mprintf("\t\t\tStep size is : %i\n",step_size);
    //mprintf("\t\t\ty is : %f\n",y[i]);
    //mprintf("\t\t\ty2 is : %f\n",y2[i]);
  }
}

/** Called once after CMAP parameters have been read. This populates
  * various partial derivatives in the cmapParameter%{dPsi,dPhi,dPsi_dPhi}
  * object. It does this using a cubic spline on the read in CMAP grid.
  *
  * Later, this information is used to evaluate an arbitary phi/phi angle
  * for a given crossterm.
  */
int CMAP::generate_cmap_derivatives(Topology const& topIn)
{
  if (!topIn.HasCmap()) {
    mprinterr("Internal Error: CMAP::generate_cmap_derivatives(): Topology '%s' does not have CMAP parameters.\n", topIn.c_str());
    return 1;
  }
  // Allocate partial derivative matrices
  cmap_dPhi_.resize( topIn.CmapGrid().size() );
  cmap_dPsi_.resize( topIn.CmapGrid().size() );
  for (unsigned int idx = 0; idx != topIn.CmapGrid().size(); idx++) {
    unsigned int res = topIn.CmapGrid()[idx].Resolution();
    cmap_dPhi_[idx].resize( res, res );
    cmap_dPsi_[idx].resize( res, res );
  } 
  // Loop over all grids
  for (CmapGridArray::const_iterator grid = topIn.CmapGrid().begin();
                                     grid != topIn.CmapGrid().end(); ++grid)
  {
    // gidx is the cmap grid index
    long int gidx = grid - topIn.CmapGrid().begin();
    int res = grid->Resolution();
    int halfRes = res / 2;
    int twoRes = res * 2;
    int step_size = 360 / res;

    std::vector<double> tmpy(twoRes);
    std::vector<double> tmpy2(twoRes);

    // 1) calculate dE/dPhi
    for (int row = 0; row < res; row++)
    {
      // Step up one row each cycle, splining across all columns

      // Fill an *extended* tmp array (tmpy) for with CMAP values
      // for the 1D splining like CHARMM.
      // It is possible CHARMM does this to avoid edge issues.

      unsigned int k=0; //interal offset counter for tmp array
      mprintf("DEBUG: Row %6i from %6i to %6i\n", row, -halfRes, res+halfRes-1);
      for (int col = -halfRes; col < res+halfRes; col++) { // -12 to 35
        tmpy[k] = grid->Grid(col, row);
        mprintf("DEBUG:\t\t%6i%12.4f\n", k, tmpy[k]);
        k++;
      }

      // Calculate spline coeffients (tmpy2) for each of the 1D
      // horizontal rows in the CMAP table
      generate_cubic_spline(twoRes, step_size, tmpy, tmpy2);
     
      // Calculate %dPhi for using each row
      // of energies and corresponding splines in tmpy and tmpy2
      for (int j = 0; j < res; j++) {
        // offset array passed to evaluate_cubic_spline
        double dPhi = evaluate_cubic_spline(step_size, tmpy, tmpy2, j+halfRes); //gbl_cmap_dPhi(i,row,j))
        cmap_dPhi_[gidx].setElement(row, j, dPhi);
      }
    } // END loop over rows
  } // END loop over grids

  return 0;
}

/** Set up CMAP-related terms */
int CMAP::Setup_CMAP_Ene(Topology const& topIn) {
  if (!topIn.HasCmap()) {
    mprintf("Warning: No CMAP parameters in '%s'\n", topIn.c_str());
    return 0;
  }
  if (generate_cmap_derivatives(topIn)) {
    mprinterr("Error: Could not set up CMAP derivatives.\n");
    return 1;
  }

  return 0;
}
