#include "CMAP.h"
#include "../CharMask.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"
#include "../ParameterTypes.h"
#include "../Topology.h"
#include "../TorsionRoutines.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
CMAP::CMAP() :
  cmapGridPtr_(0),
  selected_cmaps_(0),
  all_cmaps_(0)
{}

/** DESTRUCTOR */
CMAP::~CMAP() {
  if (selected_cmaps_ != 0) delete selected_cmaps_;
}

/** The weight matrix */
const int CMAP::wt_[16][16] = {
  {1, 0, -3,  2, 0, 0,  0,  0, -3,  0,  9, -6,  2,  0, -6,  4},
  {0, 0,  0,  0, 0, 0,  0,  0,  3,  0, -9,  6, -2,  0,  6, -4},
  {0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  9, -6,  0,  0, -6,  4},
  {0, 0,  3, -2, 0, 0,  0,  0,  0,  0, -9,  6,  0,  0,  6, -4},
  {0, 0,  0,  0, 1, 0, -3,  2, -2,  0,  6, -4,  1,  0, -3,  2},
  {0, 0,  0,  0, 0, 0,  0,  0, -1,  0,  3, -2,  1,  0, -3,  2},
  {0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  2,  0,  0,  3, -2},
  {0, 0,  0,  0, 0, 0,  3, -2,  0,  0, -6,  4,  0,  0,  3, -2},
  {0, 1, -2,  1, 0, 0,  0,  0,  0, -3,  6, -3,  0,  2, -4,  2},
  {0, 0,  0,  0, 0, 0,  0,  0,  0,  3, -6,  3,  0, -2,  4, -2},
  {0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  3,  0,  0,  2, -2},
  {0, 0, -1,  1, 0, 0,  0,  0,  0,  0,  3, -3,  0,  0, -2,  2},
  {0, 0,  0,  0, 0, 1, -2,  1,  0, -2,  4, -2,  0,  1, -2,  1},
  {0, 0,  0,  0, 0, 0,  0,  0,  0, -1,  2, -1,  0,  1, -2,  1},
  {0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  1, -1,  0,  0, -1,  1},
  {0, 0,  0,  0, 0, 0, -1,  1,  0,  0,  2, -2,  0,  0, -1,  1}
};

/// DEBUG - print gradients
static inline void print_grad(Vec3 const& dA, Vec3 const& dB, Vec3 const& dC, Vec3 const& dD)
{
  for (int i = 0; i < 3; i++)
    mprintf("DEBUG:\t\tdX %12.4f%12.4f%12.4f%12.4f\n", dA[i], dB[i], dC[i], dD[i]);
}

/** Calculate CMAP energy */
/*double CMAP::Ene_CMAP(CmapArray const& Cmaps, Frame const& frameIn)
const
{
    Vec3 dPhi_dijkl[4];
    Vec3 dPsi_djklm[4];
    double dPhi, dPsi;
  return get_cmap_energy(Cmaps, frameIn, dPhi, dPsi, dPhi_dijkl, dPsi_djklm);
}*/

/** Calculate CMAP energy and partial derivatives */
/*double CMAP::get_cmap_energy(CmapArray const& Cmaps, Frame const& frameIn,
                             double& dPhi, double& dPsi,
                             Vec3(&dPhi_dijkl)[4], Vec3(&dPsi_djklm)[4])*/

/** Calculate CMAP energy */
double CMAP::Ene_CMAP(Frame const& frameIn)
const
{
  if ( all_cmaps_ != 0)
    return Ene_CMAP( *all_cmaps_, frameIn );
  else
    return Ene_CMAP( *selected_cmaps_, frameIn );
}

/** Calculate CMAP energy */
double CMAP::Ene_CMAP(CmapArray const& Cmaps, Frame const& frameIn)
const
{
  double ene_cmap = 0.0;
  Vec3 dPhi_dijkl[4];
  Vec3 dPsi_djklm[4];
  double dPhi, dPsi;

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
    double cosphi_ijkl, sinphi_ijkl;
    Torsion_and_part_deriv( ixyz, jxyz, kxyz, lxyz,
                            dPhi_dijkl[0], dPhi_dijkl[1], dPhi_dijkl[2], dPhi_dijkl[3],
                            cosphi_ijkl, sinphi_ijkl );
    print_grad(dPhi_dijkl[0], dPhi_dijkl[1], dPhi_dijkl[2], dPhi_dijkl[3]);
//    mprintf("%30.15f\n", acos(cosphi_ijkl));
    double phi = copysign(acos(cosphi_ijkl),sinphi_ijkl) * Constants::RADDEG;
    mprintf("DEBUG: Dihedral 1 %i %i %i %i = %30.15f%30.15f%30.15f\n", cmap->A1()+1, cmap->A2()+1, cmap->A3()+1, cmap->A4()+1, phi, cosphi_ijkl, sinphi_ijkl);

    // Calculate the dihedral angle (psi) and the derivatives of the
    // four coordinates with respect to psi. Remember this subroutine is
    // operating in radians.
    double cospsi_jklm, sinpsi_jklm;
    Torsion_and_part_deriv( jxyz, kxyz, lxyz, mxyz,
                            dPsi_djklm[0], dPsi_djklm[1], dPsi_djklm[2], dPsi_djklm[3],
                            cospsi_jklm, sinpsi_jklm );
    print_grad(dPsi_djklm[0], dPsi_djklm[1], dPsi_djklm[2], dPsi_djklm[3]);
    double psi = copysign(acos(cospsi_jklm),sinpsi_jklm) * Constants::RADDEG;
    mprintf("DEBUG: Dihedral 2 %i %i %i %i = %30.15f%30.15f%30.15f\n", cmap->A2()+1, cmap->A3()+1, cmap->A4()+1, cmap->A5()+1, psi, cospsi_jklm, sinpsi_jklm);

    ene_cmap += charmm_calc_cmap_from_phi_psi(phi, psi, cmap->Idx(), dPhi, dPsi);
  } // END loop over CMAPs
  return ene_cmap;
}

/** Calculate CMAP energy and force */
double CMAP::Ene_Frc_CMAP(Frame const& frameIn)
const
{
  if ( all_cmaps_ != 0)
    return Ene_Frc_CMAP( *all_cmaps_, frameIn );
  else
    return Ene_Frc_CMAP( *selected_cmaps_, frameIn );
}

/** Calculate CMAP energy and force */
double CMAP::Ene_Frc_CMAP(CmapArray const& Cmaps, Frame const& frameIn)
const
{
  double ene_cmap = 0.0;
  Vec3 dPhi_dijkl[4];
  Vec3 dPsi_djklm[4];
  double dPhi, dPsi;

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
    double cosphi_ijkl, sinphi_ijkl;
    Torsion_and_part_deriv( ixyz, jxyz, kxyz, lxyz,
                            dPhi_dijkl[0], dPhi_dijkl[1], dPhi_dijkl[2], dPhi_dijkl[3],
                            cosphi_ijkl, sinphi_ijkl );
    print_grad(dPhi_dijkl[0], dPhi_dijkl[1], dPhi_dijkl[2], dPhi_dijkl[3]);
//    mprintf("%30.15f\n", acos(cosphi_ijkl));
    double phi = copysign(acos(cosphi_ijkl),sinphi_ijkl) * Constants::RADDEG;
    mprintf("DEBUG: Dihedral 1 %i %i %i %i = %30.15f%30.15f%30.15f\n", cmap->A1()+1, cmap->A2()+1, cmap->A3()+1, cmap->A4()+1, phi, cosphi_ijkl, sinphi_ijkl);

    // Calculate the dihedral angle (psi) and the derivatives of the
    // four coordinates with respect to psi. Remember this subroutine is
    // operating in radians.
    double cospsi_jklm, sinpsi_jklm;
    Torsion_and_part_deriv( jxyz, kxyz, lxyz, mxyz,
                            dPsi_djklm[0], dPsi_djklm[1], dPsi_djklm[2], dPsi_djklm[3],
                            cospsi_jklm, sinpsi_jklm );
    print_grad(dPsi_djklm[0], dPsi_djklm[1], dPsi_djklm[2], dPsi_djklm[3]);
    double psi = copysign(acos(cospsi_jklm),sinpsi_jklm) * Constants::RADDEG;
    mprintf("DEBUG: Dihedral 2 %i %i %i %i = %30.15f%30.15f%30.15f\n", cmap->A2()+1, cmap->A3()+1, cmap->A4()+1, cmap->A5()+1, psi, cospsi_jklm, sinpsi_jklm);

    ene_cmap += charmm_calc_cmap_from_phi_psi(phi, psi, cmap->Idx(), dPhi, dPsi);

    // Do force calc
    // Convert over to degrees per interval
    CmapGridType const& cmapGrid = (*(cmapGridPtr_))[cmap->Idx()];
    double rad_to_deg_coeff = Constants::RADDEG * 1 / (double)(360.0 / cmapGrid.Resolution());
    for (unsigned int i = 0; i < 4; i++) {
      dPhi_dijkl[i] *= rad_to_deg_coeff;
      dPsi_djklm[i] *= rad_to_deg_coeff;
      // Use chain rule to obtain the energy gradient wrt to coordinate
      dPhi_dijkl[i] *= dPhi;
      dPsi_djklm[i] *= dPsi;
    }
    mprintf("FINAL PHI\n");
    for (unsigned int i = 0; i < 4; i++)
      for (unsigned int j = 0; j < 3; j++)
        mprintf("%16.8f\n", dPhi_dijkl[i][j]);
    mprintf("FINAL PSI\n");
    for (unsigned int i = 0; i < 4; i++)
      for (unsigned int j = 0; j < 3; j++)
        mprintf("%16.8f\n", dPsi_djklm[i][j]);

    std::vector<Vec3> frc(5, Vec3(0.0));
    for (int n = 0; n < 3; n++) {
      frc[0][n] = -dPhi_dijkl[0][n];
      frc[1][n] = -dPhi_dijkl[1][n] - dPsi_djklm[0][n];
      frc[2][n] = -dPhi_dijkl[2][n] - dPsi_djklm[1][n];
      frc[3][n] = -dPhi_dijkl[3][n] - dPsi_djklm[2][n];
      frc[4][n] =                   - dPsi_djklm[3][n];
      mprintf("i %3i %16.8f\n", n, frc[0][n]);
      mprintf("j %3i %16.8f\n", n, frc[1][n]);
      mprintf("k %3i %16.8f\n", n, frc[2][n]);
      mprintf("l %3i %16.8f\n", n, frc[3][n]);
      mprintf("m %3i %16.8f\n", n, frc[4][n]);
    }
  }
  return ene_cmap;
}

/// Used to mimic Fortran modulo with REAL arguments
static inline double dmodulo(double a, double p) {
  return ( a - floor(a / p) * p );
}

static inline void print_stencil(const char* desc, const double (&matrix)[2][2])
{
  mprintf("%s:\n", desc);
  mprintf("%9.6f %9.6f\n", matrix[1][0], matrix[1][1]);
  mprintf("%9.6f %9.6f\n", matrix[0][0], matrix[0][1]);
  mprintf("\n");
}

static inline void flatten_stencil(double* out, const double(&matrix)[2][2], double step_size)
{
  out[0] = matrix[0][0] * step_size;
  out[1] = matrix[0][1] * step_size;
  out[2] = matrix[1][1] * step_size;
  out[3] = matrix[1][0] * step_size;
}

/** Calculate the CMAP energy given psi,phi and the cmap parameter */
double CMAP::charmm_calc_cmap_from_phi_psi(double phi, double psi, int cidx, double& dPhi, double& dPsi)
const
{
  CmapGridType const& cmapGrid = (*(cmapGridPtr_))[cidx];
  static const int gridOrigin = -180; ///< Where the 2D grid starts in degrees
  int step_size = 360 / cmapGrid.Resolution();
  double ene_cmap = 0.0;
  dPhi = 0.0;
  dPsi = 0.0;
  // Work out nearest complete grid point on the CMAP grid from
  // phi and psi and use this to form a 2x2 stencil
  mprintf("x= %16.8f  y= %16.8f\n", (phi - gridOrigin)/(step_size),(psi - gridOrigin)/(step_size) );
  int x = (int)( (phi - gridOrigin)/(step_size) );
  int y = (int)( (psi - gridOrigin)/(step_size) );
  mprintf("x= %6i  y= %6i\n", x, y);

  // Work out the fraction of the CMAP grid step that the interpolated
  // point takes up.
  // This will give the remainder part and then it is divided by the step size
  double phiFrac = dmodulo((phi - gridOrigin), step_size) / step_size;
  double psiFrac = dmodulo((psi - gridOrigin), step_size) / step_size;
  mprintf("phiFrac %16.8f  psiFrac %16.8f\n", phiFrac, psiFrac);

  double E_stencil[2][2];
  double dPhi_stencil[2][2];
  double dPsi_stencil[2][2];
  double dPhi_dPsi_stencil[2][2];

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      E_stencil[i][j] = cmapGrid.Grid(x+j, y+i);

      dPhi_stencil[i][j] = cmap_dPhi_[cidx].element_wrapped(y+i, x+j);

      dPsi_stencil[i][j] = cmap_dPsi_[cidx].element_wrapped(y+i, x+j);

      dPhi_dPsi_stencil[i][j] = cmap_dPhi_dPsi_[cidx].element_wrapped(y+i, x+j);
    }
  }

  print_stencil("CMAP", E_stencil);
  print_stencil("dPhi", dPhi_stencil);
  print_stencil("dPsi", dPsi_stencil);
  print_stencil("dPhi_dPsi", dPhi_dPsi_stencil);

  // Convert the 2x2 stencils into a 1D array for processing
  // by weight_stencil.
  // The array starts at the bottom left of the 2x2, working
  // around counterclockwise:
  //
  //          4 3
  //          1 2
  //
  //There may be a cleaner/better way of doing this
  //double E_stencil_1D[4];
  //double dPhi_stencil_1D[4];
  //double dPsi_stencil_1D[4];
  //double dPhi_dPsi_stencil_1D[4];

  std::vector<double> E_stencil_1D(16);
  double* dptr = &(E_stencil_1D[0]);
  flatten_stencil(dptr, E_stencil, 1);
  //mprintf("%9.6f %9.6f %9.6f %9.6f\n", E_stencil_1D[0], E_stencil_1D[1], E_stencil_1D[2], E_stencil_1D[3]);
  flatten_stencil(dptr+4, dPhi_stencil, step_size);
  flatten_stencil(dptr+8, dPsi_stencil, step_size);
  flatten_stencil(dptr+12, dPhi_dPsi_stencil, step_size*step_size);
  for (std::vector<double>::const_iterator it = E_stencil_1D.begin(); it != E_stencil_1D.end(); ++it)
    mprintf("%16.8f\n", *it);
 
/*
  E_stencil_1D(0) = E_stencil(0,0)
  E_stencil_1D(1) = E_stencil(0,1)
  E_stencil_1D(2) = E_stencil(1,1)
  E_stencil_1D(3) = E_stencil(1,0)

  !DEBUG
  !write(6,'(f9.6)'),E_stencil_1D

  dPhi_stencil_1D(1) = dPhi_stencil(1,1)
  dPhi_stencil_1D(2) = dPhi_stencil(1,2)
  dPhi_stencil_1D(3) = dPhi_stencil(2,2)
  dPhi_stencil_1D(4) = dPhi_stencil(2,1)

  dPsi_stencil_1D(1) = dPsi_stencil(1,1)
  dPsi_stencil_1D(2) = dPsi_stencil(1,2)
  dPsi_stencil_1D(3) = dPsi_stencil(2,2)
  dPsi_stencil_1D(4) = dPsi_stencil(2,1)

  dPhi_dPsi_stencil_1D(1) = dPhi_dPsi_stencil(1,1)
  dPhi_dPsi_stencil_1D(2) = dPhi_dPsi_stencil(1,2)
  dPhi_dPsi_stencil_1D(3) = dPhi_dPsi_stencil(2,2)
  dPhi_dPsi_stencil_1D(4) = dPhi_dPsi_stencil(2,1)

  !Weight d{Psi,dPhi,dPhidPsi}_stencil
  call weight_stencil(gbl_cmap_grid_step_size(cmap_idx), &
                      E_stencil_1D,              &
                      dPhi_stencil_1D,           &
                      dPsi_stencil_1D,           &
                      dPhi_dPsi_stencil_1D,      &
                      c )
*/
  // Generate coefficients for bicubic interpolation
  // NOTE: This is based on routine weight_stencil in cmap.F90
 
  std::vector<double> coeff(16);
  for (int i = 0; i < 16; i++) {
    coeff[i] = 0.0;
    for (int k = 0; k < 16; k++) {
      coeff[i] += wt_[k][i] * E_stencil_1D[k];
    }
  }
  for (std::vector<double>::const_iterator it = coeff.begin(); it != coeff.end(); ++it)
    mprintf("%16.8f\n", *it);
/*
  call bicubic_interpolation(phiFrac, psiFrac, c, E, dphi, dpsi)
*/
  // Bicubic interpolation
  // NOTE: This is based on routine bicubic_interpolation in cmap.F90
  std::vector<double>::const_reverse_iterator Cit = coeff.rbegin();
  std::vector<double>::const_reverse_iterator Dit = coeff.rbegin();
  for (int i = 3; i > -1; i--, Cit += 4, Dit++) {
    mprintf("BICUBIC %i %16.8f %16.8f %16.8f %16.8f\n", i, *(Cit), *(Cit+1), *(Cit+2), *(Cit+3));
    mprintf("BICUBIC2 %i %16.8f %16.8f %16.8f\n", i, *(Dit), *(Dit+4), *(Dit+8));
    double c4 = *(Cit);
    double c3 = *(Cit+1);
    double c2 = *(Cit+2);
    double c1 = *(Cit+3);
    double d4 = *(Dit);
    double d3 = *(Dit+4);
    double d2 = *(Dit+8);
    ene_cmap =  ene_cmap*phiFrac +( (c4*psiFrac + c3)*psiFrac + c2 )*psiFrac + c1;
    dPhi = dPhi*psiFrac +( 3.0*d4*phiFrac + 2.0*d3 )*phiFrac + d2;
    dPsi = dPsi*phiFrac +( 3.0*c4*psiFrac + 2.0*c3 )*psiFrac + c2;
  }
  mprintf("FINAL %16.8f %16.8f %16.8f\n", ene_cmap, dPhi, dPsi);

  return ene_cmap;
}

// -----------------------------------------------------------------------------
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
  cmap_dPhi_dPsi_.resize( topIn.CmapGrid().size() );
  for (unsigned int idx = 0; idx != topIn.CmapGrid().size(); idx++) {
    unsigned int res = topIn.CmapGrid()[idx].Resolution();
    cmap_dPhi_[idx].resize( res, res );
    cmap_dPsi_[idx].resize( res, res );
    cmap_dPhi_dPsi_[idx].resize( res, res );
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

      unsigned int k=0; // index for tmpy array
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

    // 2) calculate dE/dPsi
    for (int col = 0; col < res; col++)
    {
      // Step across one column each cycle, splining up each column

      unsigned int k=0; // index for tmpy array
      mprintf("DEBUG: Col %6i from %6i to %6i\n", col, -halfRes, res+halfRes-1);
      for (int row = -halfRes; row < res+halfRes; row++) { // -12 to 35
        tmpy[k] = grid->Grid(col, row);
        mprintf("DEBUG:\t\t%6i%12.4f\n", k, tmpy[k]);
        k++;
      }

      // Calculate spline coeffients (tmpy2) for each of the 1D
      // vertical columns in the CMAP table
      generate_cubic_spline(twoRes, step_size, tmpy, tmpy2);
     
      // Calculate %dPsi for using each column of energies and
      // corresponding splines in tmpy and tmpy2
      for (int j = 0; j < res; j++) {
        // offset array passed to evaluate_cubic_spline
        double dPsi = evaluate_cubic_spline(step_size, tmpy, tmpy2, j+halfRes); //gbl_cmap_dPsi(i,j,col))
        cmap_dPsi_[gidx].setElement(j, col, dPsi);
      }
    } // END loop over cols

    // 3) calculate d^2E/dPhidPsi
    // TODO Interpolate partitial derivative of psi; dE/dPhi ?
    for (int col = 0; col < res; col++)
    {
      unsigned int k=0;
      mprintf("DEBUG: D2 %6i from %6i to %6i\n", col, -halfRes, res+halfRes-1);
      for (int row = -halfRes; row < res+halfRes; row++) { // -12 to 35
        tmpy[k] = cmap_dPhi_[gidx].element_wrapped(row, col);
        mprintf("DEBUG:\t\t%6i%12.4f\n", k, tmpy[k]);
        k++;
      }

      generate_cubic_spline(twoRes, step_size, tmpy, tmpy2);
     
      for (int j = 0; j < res; j++) {
        double dPhiPsi = evaluate_cubic_spline(step_size, tmpy, tmpy2, j+halfRes); //gbl_cmap_dPhi_dPsi(i,j,col))
        cmap_dPhi_dPsi_[gidx].setElement(j, col, dPhiPsi);
      }
    } // END loop over cols

  } // END loop over grids

  return 0;
}

/** Set up CMAP-related terms */
int CMAP::Setup_CMAP_Ene(Topology const& topIn, CharMask const& cmask) {
  if (selected_cmaps_ != 0) delete selected_cmaps_;
  selected_cmaps_ = 0;
  all_cmaps_ = 0;

  if (!topIn.HasCmap()) {
    mprintf("Warning: No CMAP parameters in '%s'\n", topIn.c_str());
    return 0;
  }
  cmapGridPtr_ = &(topIn.CmapGrid());
  if (generate_cmap_derivatives(topIn)) {
    mprinterr("Error: Could not set up CMAP derivatives.\n");
    return 1;
  }

  CmapArray const* cmaparray = 0;
  if (cmask.All()) {
    mprintf("\tAll CMAPs will be selected.\n");
    all_cmaps_ = &(topIn.Cmap());
    cmaparray = all_cmaps_;
    mprintf("\t%zu CMAPs.\n", all_cmaps_->size());
  } else {
    selected_cmaps_ = new CmapArray();
    cmaparray = (CmapArray const*)selected_cmaps_;
    for (CmapArray::const_iterator cmap = topIn.Cmap().begin();
                                   cmap != topIn.Cmap().end(); ++cmap)
    {
      if ( cmask.AtomInCharMask(cmap->A1()) &&
           cmask.AtomInCharMask(cmap->A2()) &&
           cmask.AtomInCharMask(cmap->A3()) &&
           cmask.AtomInCharMask(cmap->A4()) &&
           cmask.AtomInCharMask(cmap->A5()) )
      {
        selected_cmaps_->push_back( *cmap );
      }
    }
    mprintf("\t%zu CMAPs.\n", selected_cmaps_->size());
  }
  for (CmapArray::const_iterator cmap = cmaparray->begin(); cmap != cmaparray->end(); ++cmap)
  {
    mprintf("\t\tSelecting CMAP %s - %s - %s - %s - %s\n",
            topIn.AtomMaskName(cmap->A1()).c_str(),
            topIn.AtomMaskName(cmap->A2()).c_str(),
            topIn.AtomMaskName(cmap->A3()).c_str(),
            topIn.AtomMaskName(cmap->A4()).c_str(),
            topIn.AtomMaskName(cmap->A5()).c_str());
  }
  return 0;
}
