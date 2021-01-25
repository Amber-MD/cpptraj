#ifndef INC_ENERGYKERNEL_NONBOND_SIMPLE_H
#define INC_ENERGYKERNEL_NONBOND_SIMPLE_H
#include <cmath> //sqrt
#include "Frame.h"
#include "CharMask.h"
/// Simple LJ and electrostatics
template <class REAL> class EnergyKernel_NonBond_Simple {
  public:
     /// CONSTRUCTOR
     EnergyKernel_NonBond_Simple() {}
     /// Energy/forces
     void Calc_F_E(Frame&, int, int, double, double, double, CharMask const&, double&, double&);
};

template<class REAL> 
void EnergyKernel_NonBond_Simple<REAL>::Calc_F_E(Frame& frameIn, int idx, int jdx,
                                           double LJA, double LJB, double qiqj,
                                           CharMask const& maskIn,
                                           double& E_vdw, double& E_elec)
{
  const double* XYZ0 = frameIn.XYZ( idx );
  const double* XYZ1 = frameIn.XYZ( jdx );
  REAL rx = XYZ0[0] - XYZ1[0];
  REAL ry = XYZ0[1] - XYZ1[1];
  REAL rz = XYZ0[2] - XYZ1[2];
  REAL rij2 = rx*rx + ry*ry + rz*rz;
  if (rij2 > 0) {
    REAL rij = sqrt( rij2 );
    // VDW
    REAL r2    = 1.0 / rij2;
    REAL r6    = r2 * r2 * r2;
    REAL r12   = r6 * r6;
    REAL f12   = LJA * r12;  // A/r^12
    REAL f6    = LJB * r6;   // B/r^6
    REAL e_vdw = f12 - f6;   // (A/r^12)-(B/r^6)
    //mprintf("DBG:\t\t%8i %8i %12.4f\n", e_vdw);
    E_vdw += e_vdw;
    // VDW force
    REAL fvdw = ((12*f12) - (6*f6)) * r2; // (12A/r^13)-(6B/r^7)
    REAL dfx = rx * fvdw;
    REAL dfy = ry * fvdw;
    REAL dfz = rz * fvdw;
    // COULOMB
    REAL e_elec = 1.0 * (qiqj / rij); // 1.0 is electrostatic constant, not really needed
    E_elec += e_elec;
    // COULOMB force
    REAL felec = e_elec / rij; // kes * (qiqj / r) * (1/r)
    dfx += rx * felec;
    dfy += ry * felec;
    dfz += rz * felec;
    // Apply forces
    if (maskIn.AtomInCharMask(idx)) {
      double* fxyz = frameIn.fAddress() + (3*idx);
      fxyz[0] += dfx;
      fxyz[1] += dfy;
      fxyz[2] += dfz;
    }
    if (maskIn.AtomInCharMask(jdx)) {
      double* fxyz = frameIn.fAddress() + (3*jdx);
      fxyz[0] -= dfx;
      fxyz[1] -= dfy;
      fxyz[2] -= dfz;
    }
  } // END rij > 0
}
#endif
