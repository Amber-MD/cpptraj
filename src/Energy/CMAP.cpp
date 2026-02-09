#include "CMAP.h"
#include "../CpptrajStdio.h" // DEBUG
#include "../Frame.h"
#include "../ParameterTypes.h"
#include "../TorsionRoutines.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
CMAP::CMAP() {}

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
    double psi = copysign(acos(cospsi_jklm),sinpsi_jklm) * Constants::RADDEG;
    mprintf("DEBUG: Dihedral 2 %i %i %i %i = %g deg.\n", cmap->A2()+1, cmap->A3()+1, cmap->A4()+1, cmap->A5()+1, psi);
  }
  return ene_cmap;
}
