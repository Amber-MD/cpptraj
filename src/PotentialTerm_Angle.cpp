#include "PotentialTerm_Angle.h"
#include "Topology.h"
#include "CharMask.h"
#include "EnergyArray.h"
//#incl ude "CpptrajStdio.h" 

void PotentialTerm_Angle::addAngles(AngleArray const& angles, CharMask const& maskIn)
{
  for (AngleArray::const_iterator ang = angles.begin(); ang != angles.end(); ++ang)
  {
    if (maskIn.AtomInCharMask( ang->A1() ) ||
        maskIn.AtomInCharMask( ang->A2() ) ||
        maskIn.AtomInCharMask( ang->A3() ))
    {
      //mprintf("DEBUG: Angle %i to %i to %i\n", ang->A1()+1, ang->A2()+1, ang->A3()+1);
      activeAngs_.push_back( *ang );
    }
  }
}

/** Set up Hooke's law angle term. */
int PotentialTerm_Angle::SetupTerm(Topology const& topIn, Box const& boxIn,
                                   CharMask const& maskIn, EnergyArray& Earray)
{
  activeAngs_.clear();
  addAngles( topIn.Angles(),  maskIn );
  addAngles( topIn.AnglesH(), maskIn );

  angParm_ = &(topIn.AngleParm());
  Eang_ = Earray.AddType( EnergyArray::E_ANGLE );

  return 0;
}

/** Calculate Hooke's law angle force.
  * NOTE: Code adapted from AmberTools 20 SFF, sff2.c eangl2()
  */
void PotentialTerm_Angle::CalcForce(Frame& frameIn, CharMask const& maskIn) const {
  *Eang_ = 0.0;

  // TODO rename these, fix indices
  double s[7];
  double dc[7];
  double dr[10];
  //mprintf("FCALC\n");
  for (AngleArray::const_iterator ang = activeAngs_.begin(); ang != activeAngs_.end(); ++ang)
  {
    AngleParmType AP = (*angParm_)[ ang->Idx() ];
    const double* XYZ1 = frameIn.XYZ( ang->A1() );
    const double* XYZ2 = frameIn.XYZ( ang->A2() );
    const double* XYZ3 = frameIn.XYZ( ang->A3() );
/*
    // TODO just use CalcAngle?
    Vec3 v2_1 = Vec3( XYZ1[0] - XYZ2[0],
                      XYZ1[1] - XYZ2[1],
                      XYZ1[2] - XYZ2[2] );

    Vec3 v2_3 = Vec3( XYZ3[0] - XYZ2[0],
                      XYZ3[1] - XYZ2[1],
                      XYZ3[2] - XYZ2[2] );

    Vec3 v1_3 = Vec3( XYZ3[0] - XYZ1[0],
                      XYZ3[1] - XYZ1[1],
                      XYZ3[2] - XYZ1[2] );

    double mag_v2_1 = v2_1.Magnitude2();
    double ang_in_rad = 0;
    if (mag_v2_1 > Constants::SMALL) {
      double mag_v2_3 = v2_3.Magnitude2();
      if (mag_v2_3 > Constants::SMALL) {
        double cos_ang = (v2_1 * v2_3) / sqrt(mag_v2_1 * mag_v2_3);
        if (cos_ang > 1.0)
          cos_ang = 1.0;
        else if (cos_ang < -1.0)
          cos_ang = -1.0;
        ang_in_rad = acos( cos_ang );
      }
    }
    mprintf("FCALC ANGLE %20.10f\n", ang_in_rad);
*/

    double x1 = XYZ1[0];
    double y1 = XYZ1[1];
    double z1 = XYZ1[2];

    double x2 = XYZ2[0];
    double y2 = XYZ2[1];
    double z2 = XYZ2[2];

    double x3 = XYZ3[0];
    double y3 = XYZ3[1];
    double z3 = XYZ3[2];

    s[1] = x1 - x2;
    s[2] = y1 - y2;
    s[3] = z1 - z2;
    s[4] = x3 - x2;
    s[5] = y3 - y2;
    s[6] = z3 - z2;
    double rij2 = s[1] * s[1] + s[2] * s[2] + s[3] * s[3];
    double rkj2 = s[4] * s[4] + s[5] * s[5] + s[6] * s[6];
    double rij = sqrt(rij2);
    double rkj = sqrt(rkj2);
    double rrik = rij * rkj;
    double c = (s[1] * s[4] + s[2] * s[5] + s[3] * s[6]) / rrik;

    c = c > 1.0 ? 1.0 : c;
    c = c < -1.0 ? -1.0 : c;
    double theta = acos(c);
    //mprintf("FCALC ANGLE %20.10f\n", theta);
    double dtheta = theta - AP.Teq();

    // df and ddf are derivatives of E with respect to c == costheta:
    double df = dtheta * AP.Tk();
    //*double ddf = 2. * AP.Tk();*
    double e = df * dtheta;
    *Eang_ += e;

    // Calculate the sine then limit small values of the sine to the range
    // -10**-3..+10**-3 so as to avoid division by zero.
    double st = sin(theta);
    if (st > 0 && st < 1.e-3)
       st = 1.e-3;
    else if (st < 0 && st > -1.e-3)
       st = -1.e-3;
    df *= 2.0;                // check this!!  

    // dc = derivative of c with respect to cartesian differences:
    for (int i = 1; i <= 3; i++) {
       dc[i] = (s[i + 3] / rkj - c * s[i] / rij) / rij;
       dc[i + 3] = (s[i] / rij - c * s[i + 3] / rkj) / rkj;
    }
 /*
    // get ddc = second derivative of c with respect to 
    // cartesian differences:
    difang(c, s, dc, ddc);

    // change ddc to second. der. of theta with respect to 
    // cartesian differences:
    st2c = c / (st * st);
    for (i = 1; i <= 6; i++) {
       for (j = i; j <= 6; j++) {
          ddc[i][j] = -(ddc[i][j] + dc[i] * dc[j] * st2c) / st;
          ddc[j][i] = ddc[i][j];
       }
    }
 */

    // change dc to derivative of theta w/ respect to cartesian
    // differences:
    for (int i = 1; i <= 6; i++) {
       dc[i] = -dc[i] / st;
    }

    // dr will hold -derivates of theta w/ respect to cartesians:
    for (int i = 1; i <= 3; i++) {
       dr[i]     = -dc[i];
       dr[i + 6] = -dc[i + 3];
       dr[i + 3] =  dc[i] + dc[i + 3];
    }

    // update the forces:  
    if (maskIn.AtomInCharMask(ang->A1())) {
      double* frc = frameIn.fAddress() + ang->A1()*3;
      frc[0] += df * dr[1];
      frc[1] += df * dr[2];
      frc[2] += df * dr[3];
    }

    if (maskIn.AtomInCharMask(ang->A2())) {
      double* frc = frameIn.fAddress() + ang->A2()*3;
      frc[0] += df * dr[4];
      frc[1] += df * dr[5];
      frc[2] += df * dr[6];
    }

    if (maskIn.AtomInCharMask(ang->A3())) {
      double* frc = frameIn.fAddress() + ang->A3()*3;
      frc[0] += df * dr[7];
      frc[1] += df * dr[8];
      frc[2] += df * dr[9];
    }
  } // END loop over angles
}
