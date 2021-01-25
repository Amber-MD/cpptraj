#include "PotentialTerm_Dihedral.h"
#include "Topology.h"
#include "CharMask.h"
#include "EnergyArray.h"
#include "CpptrajStdio.h" 

void PotentialTerm_Dihedral::addDihedrals(DihedralArray const& dihedrals, CharMask const& maskIn)
{
  for (DihedralArray::const_iterator dih = dihedrals.begin(); dih != dihedrals.end(); ++dih)
  {
    if (maskIn.AtomInCharMask( dih->A1() ) ||
        maskIn.AtomInCharMask( dih->A2() ) ||
        maskIn.AtomInCharMask( dih->A3() ) ||
        maskIn.AtomInCharMask( dih->A4() ))
    {
      mprintf("DEBUG: Dihedral %i to %i to %i to %i\n", dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1);
      activeDihs_.push_back( *dih );
    }
  }
}

/** Set up truncated Fourier series dihedral term. */
int PotentialTerm_Dihedral::SetupTerm(Topology const& topIn, Box const& boxIn,
                                      CharMask const& maskIn, EnergyArray& Earray)
{
  activeDihs_.clear();
  addDihedrals( topIn.Dihedrals(),  maskIn );
  addDihedrals( topIn.DihedralsH(), maskIn );

  dihParm_ = &(topIn.DihedralParm());
  Edih_ = Earray.AddType( EnergyArray::E_DIHEDRAL );

  return 0;
}

/** Calculate truncated Fourier series dihedral term. 
  * NOTE: Code adapted from AmberTools 20 SFF, sff2.c ephi2()
  */
void PotentialTerm_Dihedral::CalcForce(Frame& frameIn, CharMask const& maskIn) const {
  *Edih_ = 0.0;

  for (DihedralArray::const_iterator dih = activeDihs_.begin(); dih != activeDihs_.end(); ++dih)
  {
    DihedralParmType DP = (*dihParm_)[ dih->Idx() ];
    const double* XYZ1 = frameIn.XYZ( dih->A1() );
    const double* XYZ2 = frameIn.XYZ( dih->A2() );
    const double* XYZ3 = frameIn.XYZ( dih->A3() );
    const double* XYZ4 = frameIn.XYZ( dih->A4() );

    double xij = XYZ1[0] - XYZ2[0];
    double yij = XYZ1[1] - XYZ2[1];
    double zij = XYZ1[2] - XYZ2[2];

    double xkj = XYZ3[0] - XYZ2[0];
    double ykj = XYZ3[1] - XYZ2[1];
    double zkj = XYZ3[2] - XYZ2[2];

    double xkl = XYZ3[0] - XYZ4[0];
    double ykl = XYZ3[1] - XYZ4[1];
    double zkl = XYZ3[2] - XYZ4[2];

    double dx = yij * zkj - zij * ykj;
    double dy = zij * xkj - xij * zkj;
    double dz = xij * ykj - yij * xkj;

    double gx = zkj * ykl - ykj * zkl;
    double gy = xkj * zkl - zkj * xkl;
    double gz = ykj * xkl - xkj * ykl;

    double bi = dx * dx + dy * dy + dz * dz;
    double bk = gx * gx + gy * gy + gz * gz;
    double ct = dx * gx + dy * gy + dz * gz;

    //    ----- approximate if linear dihedral 
    //    ----- assumes force constant is eventually zero -----

    if ( bi <= 0.01 ) bi = 0.01;
    if ( bk <= 0.01 ) bk = 0.01;

    bi = sqrt(bi);
    bk = sqrt(bk);
    double z1 = 1. / bi;
    double z2 = 1. / bk;
    ct = ct * z1 * z2;
    ct = ct > 1.0 ? 1.0 : ct;
    ct = ct < -1.0 ? -1.0 : ct;

    /*     ----- ct value above is actually -cosphi; here we change
     *           its sign -----
     */
    ct = -ct;

    /*     ----- calculate the energy and derivatives -----


     *     ----- get df = first der. of potential w/respect to cosphi; and
     *           ddf = second der. of potntial w/respect to cosphi -----

     *           the torsional potential is assumed to have the form:
     *            e = pk(ic) * (1.0+phase*cos(pn(ic)*phi)
     *            where phase = 1.0 or -1.0, and pn = 1,2,3,4, or 6

     *     ----- energy terms for dihedrals are expressed in terms of
     *           cosphi in order to eliminate problems for planar angles ----
     */

    // multi-term // TODO use Constants PI?
    double df;
    if ((fabs(DP.Phase() - 3.142) > 0.01) && (fabs(DP.Phase()) > 0.01)) {
      /*   here we have special code for phases that are not zero or pi */

      double phi = acos(ct);

      /*
         now calculate sin(phi) because cos(phi) is symmetric, so
         we can decide between +-phi.
       */

      double ux = -yij * zkj + zij * ykj;
      double uy = -zij * xkj + xij * zkj;
      double uz = -xij * ykj + yij * xkj;

      double vx = -ykj * zkl + zkj * ykl;
      double vy = -zkj * xkl + xkj * zkl;
      double vz = -xkj * ykl + ykj * xkl;

      double dx1 = uy * vz - uz * vy;
      double dy1 = uz * vx - ux * vz;
      double dz1 = ux * vy - uy * vx;

      dx1 = dx1*xkj + dy1*ykj + dz1*zkj;
      if (dx1 < 0.0)
         phi = -phi;

      double delta = DP.Pn()*phi - DP.Phase();
      double e = DP.Pk()*(1.0 + cos(delta));
      *Edih_ += e;
      double yy = sin(phi);

      if (fabs(yy) > 0.001) {
        df = DP.Pk()*DP.Pn()*sin(delta)/yy;
        //ddf = -Pk[atyp]*Pn[atyp]*( Pn[atyp]*cos(delta) -
        //        ct*sin(delta)/yy )/(yy*yy);
      } else {
        /* if sin(phi) happens to be at zero or pi, adjust it slightly */
        if (phi > -1. && phi < 1.) {  /* set sin(phi) = 0.001 */
          df = DP.Pk()*DP.Pn()*sin(delta) * 1000.;
          //ddf = -Pk[atyp]*Pn[atyp]*1000000. * ( Pn[atyp]*cos(delta) -
          //        ct*1000.*sin(delta) );
        } else {  /* set sin(phi) = -0.001  */
          df = -DP.Pk()*DP.Pn()*sin(delta) * 1000.;
          //ddf = -Pk[atyp]*Pn[atyp]*1000000. * ( Pn[atyp]*cos(delta) +
          //        ct*1000.*sin(delta) );
        }
      }

    } else {   /* here is "usual" case, where the phase is 0 or pi */
      *Edih_ = 0; // FIXME

    }
  } // END loop over dihedrals
}
