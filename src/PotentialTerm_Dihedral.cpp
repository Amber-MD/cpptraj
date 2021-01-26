#include "PotentialTerm_Dihedral.h"
#include "Topology.h"
#include "CharMask.h"
#include "EnergyArray.h"
#include "EnergyKernel_NonBond_Simple.h"
#include "Constants.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
PotentialTerm_Dihedral::PotentialTerm_Dihedral() :
  PotentialTerm(DIHEDRAL),
  dihParm_(0),
  nonbond_(0),
  atoms_(0),
  Edih_(0),
  Enb14_(0),
  Eq14_(0)
{}


void PotentialTerm_Dihedral::addDihedrals(DihedralArray const& dihedrals, CharMask const& maskIn)
{
  for (DihedralArray::const_iterator dih = dihedrals.begin(); dih != dihedrals.end(); ++dih)
  {
    if (maskIn.AtomInCharMask( dih->A1() ) ||
        maskIn.AtomInCharMask( dih->A2() ) ||
        maskIn.AtomInCharMask( dih->A3() ) ||
        maskIn.AtomInCharMask( dih->A4() ))
    {
      //mprintf("DEBUG: Dihedral %i to %i to %i to %i\n", dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1);
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
  nonbond_ = &(topIn.Nonbond());
  atoms_ = &(topIn.Atoms());
  Edih_ = Earray.AddType( EnergyArray::E_DIHEDRAL );
  Enb14_ = Earray.AddType( EnergyArray::E_VDW14 );
  Eq14_ = Earray.AddType( EnergyArray::E_Q14 );

  return 0;
}

/** Calculate truncated Fourier series dihedral term. 
  * NOTE: Code adapted from AmberTools 20 SFF, sff2.c ephi2()
  */
void PotentialTerm_Dihedral::CalcForce(Frame& frameIn, CharMask const& maskIn) const {
  *Edih_ = 0.0;
  *Enb14_ = 0.0;
  *Eq14_ = 0.0;

  // TODO fix indices and names
  double t[7];
  double dc[7];
  double dr[13];
  double dtx[7][13];

  double QFAC = Constants::ELECTOAMBER * Constants::ELECTOAMBER;

  EnergyKernel_NonBond_Simple<double> NB14;
  //mprintf("FCALC\n");
  for (DihedralArray::const_iterator dih = activeDihs_.begin(); dih != activeDihs_.end(); ++dih)
  {
    DihedralParmType DP = (*dihParm_)[ dih->Idx() ];

    // Do the 1-4 terms
    if (!dih->Skip14()) {
      Atom const& A1 = (*atoms_)[dih->A1()];
      Atom const& A4 = (*atoms_)[dih->A4()];
      int nbindex = nonbond_->GetLJindex( A1.TypeIndex(), A4.TypeIndex() );
      NonbondType const& LJ = nonbond_->NBarray( nbindex );
      NB14.Calc_F_E( frameIn, dih->A1(), dih->A4(), LJ.A(), LJ.B(),
                     QFAC, A1.Charge(), A4.Charge(),
                     1.0/DP.SCNB(), 1.0/DP.SCEE(), maskIn,
                     *Enb14_, *Eq14_);
    }

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
    //mprintf("FCALC CT %20.10f\n", ct);
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

    //     ----- ct value above is actually -cosphi; here we change
    //           its sign -----
    ct = -ct;

    //     ----- calculate the energy and derivatives -----


    //     ----- get df = first der. of potential w/respect to cosphi; and
    //           ddf = second der. of potntial w/respect to cosphi -----

    //           the torsional potential is assumed to have the form:
    //            e = pk(ic) * (1.0+phase*cos(pn(ic)*phi)
    //            where phase = 1.0 or -1.0, and pn = 1,2,3,4, or 6

    //     ----- energy terms for dihedrals are expressed in terms of
    //           cosphi in order to eliminate problems for planar angles ----

    // multi-term // TODO use Constants PI?
    double df;
    if ((fabs(DP.Phase() - 3.142) > 0.01) && (fabs(DP.Phase()) > 0.01)) {
      //   here we have special code for phases that are not zero or pi

      double phi = acos(ct);

      // now calculate sin(phi) because cos(phi) is symmetric, so
      // we can decide between +-phi.

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
        // if sin(phi) happens to be at zero or pi, adjust it slightly
        if (phi > -1. && phi < 1.) {
          // set sin(phi) = 0.001
          df = DP.Pk()*DP.Pn()*sin(delta) * 1000.;
          //ddf = -Pk[atyp]*Pn[atyp]*1000000. * ( Pn[atyp]*cos(delta) -
          //        ct*1000.*sin(delta) );
        } else {
          // set sin(phi) = -0.001
          df = -DP.Pk()*DP.Pn()*sin(delta) * 1000.;
          //ddf = -Pk[atyp]*Pn[atyp]*1000000. * ( Pn[atyp]*cos(delta) +
          //        ct*1000.*sin(delta) );
        }
      }

    } else {
      // here is "usual" case, where the phase is 0 or pi
      int iper = (int)fabs(DP.Pn());
      // TODO check that Pn is not a float?
      //int iper = (int) fabs(DP.Pn()) + 0.0001;
      //assert(iper != (fabs(DP.Pn()) - 0.0001));
      //assert(iper >= 1 && iper <= 6);
      double e = 0;
      double ct2 = 0;
      switch (iper) {
      case 1:
         e = ct;
         df = 1.0;
         //ddf = 0.0;
         break;

      case 2:
         e = 2.0*ct*ct - 1.0;
         df = 4.0*ct;
         //ddf = 4.0;
         break;

      case 3:
         ct2 = ct*ct;
         e = ct*(4.0*ct2 - 3.0);
         df = 12.0*ct2 - 3.0;
         //ddf = 24.0*ct;
         break;

      case 4:
         ct2 = ct*ct;
         e = 1.0 + ct2*8.0*(ct2 - 1.0);
         df = 32.0*ct2*ct - 16.*ct;
         //ddf = 96.*ct2 -16.;
         break;

      case 5:
         ct2 = ct*ct;
         e = 16.*ct2*ct2*ct - 20.*ct2*ct + 5.*ct;
         df = 80.*ct2*ct2 - 60.*ct2 + 5.;
         //ddf = 320.*ct2*ct - 120.*ct;
         break;

      case 6:
         ct2 = ct * ct;
         e = ct2 * (ct2 * (ct2 * 32.0 - 48.0) + 18.0) - 1.0;
         df = ct * (ct2 * (ct2 * 192.0 - 192.0) + 36.0);
         //ddf = ct2 * (ct2 * 960.0 - 576.0) + 36.0;
         break;

      default:
         mprinterr("Error: PotentialTerm_Dihedral: bad periodicity: %i\n", iper);
         return;
      } // END switch over iper

      double arg;
      if (fabs(DP.Phase() - 3.142) < 0.01)
         arg = -DP.Pk();
      else
         arg = DP.Pk();

      e = DP.Pk() + arg * e;
      df = df * arg;
      //ddf = ddf * arg;
      *Edih_ += e;
    }

    t[1] = dx;
    t[2] = dy;
    t[3] = dz;

    t[4] = -gx;
    t[5] = -gy;
    t[6] = -gz;

    /*
     *     ----- now, set up array dc = first der. of cosphi w/respect
     *           to the cartesian differences t -----
     */

    double z11 = z1 * z1;
    double z12 = z1 * z2;
    double z22 = z2 * z2;
    for (int i = 1; i <= 3; i++) {
       dc[i] = t[i + 3] * z12 - ct * t[i] * z11;
       dc[i + 3] = t[i] * z12 - ct * t[i + 3] * z22;
    }

    /*
     *     ----- subroutine difang will now create array ddc which is second
     *           derivative of cosphi with respect to the t s -----
     */

    //difang(ct, t, dc, ddc);

    /* ---now set up array s, given on page 118 of cff book-- */

    double s1 = xij;
    double s2 = yij;
    double s3 = zij;
    double s4 = xkj;
    double s5 = ykj;
    double s6 = zkj;
    double s7 = -xkj;
    double s8 = -ykj;
    double s9 = -zkj;
    double s10 = -xkl;
    double s11 = -ykl;
    double s12 = -zkl;

    /*
     *     ----- set up dtx[i][j] = derivative of t[i] w/respect to x[j]
     *           see p. 120 of cff book -----
     */

    for (int i = 1; i <= 6; i++) {
       for (int j = 1; j <= 12; j++) {
          dtx[i][j] = 0.0;
       }
    }

    dtx[1][2] = s6;
    dtx[1][3] = -s5;
    dtx[1][5] = s3 - s6;
    dtx[1][6] = s5 - s2;
    dtx[1][8] = -s3;
    dtx[1][9] = s2;
    dtx[2][1] = -s6;
    dtx[2][3] = s4;
    dtx[2][4] = s6 - s3;
    dtx[2][6] = s1 - s4;
    dtx[2][7] = s3;
    dtx[2][9] = -s1;
    dtx[3][1] = s5;
    dtx[3][2] = -s4;
    dtx[3][4] = s2 - s5;
    dtx[3][5] = s4 - s1;
    dtx[3][7] = -s2;
    dtx[3][8] = s1;
    dtx[4][5] = s12;
    dtx[4][6] = -s11;
    dtx[4][8] = s9 - s12;
    dtx[4][9] = s11 - s8;
    dtx[4][11] = -s9;
    dtx[4][12] = s8;
    dtx[5][4] = -s12;
    dtx[5][6] = s10;
    dtx[5][7] = s12 - s9;
    dtx[5][9] = s7 - s10;
    dtx[5][10] = s9;
    dtx[5][12] = -s7;
    dtx[6][4] = s11;
    dtx[6][5] = -s10;
    dtx[6][7] = s8 - s11;
    dtx[6][8] = s10 - s7;
    dtx[6][10] = -s8;
    dtx[6][11] = s7;

    /*
     *     ----- set up dr array, containing -first derivative of cosphi with
     *           respect to cartesians -----
     */

    for (int i = 1; i <= 12; i++) {
       double dum = 0.0;
       for (int j = 1; j <= 6; j++) {
          dum += dc[j] * dtx[j][i];
       }
       dr[i] = -dum;
    }

    /* update the forces:  */
/*
    f[at1] -= df * dr[1];
    f[at1 + 1] -= df * dr[2];
    f[at1 + 2] -= df * dr[3];

    f[at2] -= df * dr[4];
    f[at2 + 1] -= df * dr[5];
    f[at2 + 2] -= df * dr[6];

    f[at3] -= df * dr[7];
    f[at3 + 1] -= df * dr[8];
    f[at3 + 2] -= df * dr[9];

    f[at4] -= df * dr[10];
    f[at4 + 1] -= df * dr[11];
    f[at4 + 2] -= df * dr[12];*/
    //mprintf("FCALC forces ");
    if (maskIn.AtomInCharMask(dih->A1())) {
      double* frc = frameIn.fAddress() + dih->A1()*3;
      frc[0] += df * dr[1];
      frc[1] += df * dr[2];
      frc[2] += df * dr[3];
      //mprintf("%16.8f%16.8f%16.8f",df*dr[1],df*dr[2],df*dr[3]);
    }

    if (maskIn.AtomInCharMask(dih->A2())) {
      double* frc = frameIn.fAddress() + dih->A2()*3;
      frc[0] += df * dr[4];
      frc[1] += df * dr[5];
      frc[2] += df * dr[6];
      //mprintf("%16.8f%16.8f%16.8f",df*dr[4],df*dr[5],df*dr[6]);
    }

    if (maskIn.AtomInCharMask(dih->A3())) {
      double* frc = frameIn.fAddress() + dih->A3()*3;
      frc[0] += df * dr[7];
      frc[1] += df * dr[8];
      frc[2] += df * dr[9];
      //mprintf("%16.8f%16.8f%16.8f",df*dr[7],df*dr[8],df*dr[9]);
    }

    if (maskIn.AtomInCharMask(dih->A4())) {
      double* frc = frameIn.fAddress() + dih->A4()*3;
      frc[0] += df * dr[10];
      frc[1] += df * dr[11];
      frc[2] += df * dr[12];
      //mprintf("%16.8f%16.8f%16.8f",df*dr[10],df*dr[11],df*dr[12]);
    }
    //mprintf("\n");

  } // END loop over dihedrals
}
