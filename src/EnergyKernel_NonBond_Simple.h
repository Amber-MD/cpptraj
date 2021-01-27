#ifndef INC_ENERGYKERNEL_NONBOND_SIMPLE_H
#define INC_ENERGYKERNEL_NONBOND_SIMPLE_H
#include <cmath> //sqrt
#include "Frame.h"
#include "CharMask.h"
//#incl ude "CpptrajStdio.h" // DEBUG
/// Simple LJ and electrostatics
template <class REAL> class EnergyKernel_NonBond_Simple {
  public:
     /// CONSTRUCTOR - Huge cutoff
     EnergyKernel_NonBond_Simple() : cutoff2_(999999.0) {}
     /// CONSTRUCTOR - specified cutoff^2
     EnergyKernel_NonBond_Simple(REAL cutIn) : cutoff2_(cutIn) {}
     /// Energy/forces
     void Calc_F_E(Frame&, int, int, double, double, double, double, double, double, double, CharMask const&, double&, double&);
  private:
    REAL cutoff2_; ///< Only calculate interactions within this cutoff (squared)
};

template<class REAL> 
void EnergyKernel_NonBond_Simple<REAL>::Calc_F_E(Frame& frameIn, int idx, int jdx,
                                           double LJA, double LJB,
                                           double QFAC, double qi, double qj,
                                           double enbfac, double eelfac,
                                           CharMask const& maskIn,
                                           double& E_vdw, double& E_elec)
{
  const double* XYZ0 = frameIn.XYZ( idx );
  const double* XYZ1 = frameIn.XYZ( jdx );
  REAL rx = XYZ0[0] - XYZ1[0];
  REAL ry = XYZ0[1] - XYZ1[1];
  REAL rz = XYZ0[2] - XYZ1[2];
  REAL rij2 = rx*rx + ry*ry + rz*rz;
  if (rij2 < cutoff2_) {
    REAL rij = sqrt( rij2 );
    // VDW
    REAL r2    = 1.0 / rij2;
    REAL r6    = r2 * r2 * r2;
    REAL r12   = r6 * r6;
    REAL f12   = LJA * r12;  // A/r^12
    REAL f6    = LJB * r6;   // B/r^6
    REAL e_vdw = f12 - f6;   // (A/r^12)-(B/r^6)
    //mprintf("DBG:\t\t%8i %8i %12.4f\n", e_vdw);
    E_vdw += (e_vdw * enbfac);
    // VDW force
//    REAL fvdw = ((12*f12) - (6*f6)) * enbfac * r2; // (12A/r^13)-(6B/r^7)
//    REAL dfx = rx * fvdw;
//    REAL dfy = ry * fvdw;
//    REAL dfz = rz * fvdw;
    // COULOMB
    REAL e_elec = QFAC * ((eelfac*qi*qj) / rij); // 1.0 is electrostatic constant, not really needed
    E_elec += e_elec;
    // COULOMB force
//    REAL felec = e_elec / rij; // kes * (qiqj / r) * (1/r)
//    dfx += rx * felec;
//    dfy += ry * felec;
//    dfz += rz * felec;
    // COMBINED Vdw and Coulomb force 
    REAL dfn = ( (((12*f12) - (6*f6)) * enbfac) + e_elec ) * r2;
    REAL dfx = rx * dfn;
    REAL dfy = ry * dfn;
    REAL dfz = rz * dfn;
    //mprintf("FCALC 14 %6i%6i%16.8f%16.8f%16.8f%16.8f\n", 3*idx,3*jdx,dfx, dfy, dfz,dfn);
    //mprintf("FCALC 14 %6i%6i%16.8f%16.8f%16.8f%16.8f\n", 3*idx,3*jdx, rx, ry, rz, dfn);
    //mprintf("FCALC 14 %6i%6i%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", 3*idx,3*jdx, f12, f6, enbfac, qi*qj*QFAC, 1/rij, eelfac, r2);
    //mprintf("FCALC 14 %6i%6i%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", 3*idx,3*jdx, f12, f6, enbfac, qi*qj*QFAC/rij, eelfac, r2);
    //mprintf("FCALC 14 %6i%6i%16.8f%16.8f%16.8f%16.8f%16.8f\n", 3*idx,3*jdx, f12, f6, enbfac, e_elec, r2);
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
  } // END rij < cutoff 
}
#endif
