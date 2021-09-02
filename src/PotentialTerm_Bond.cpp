#include "PotentialTerm_Bond.h"
#include "Topology.h"
#include "CharMask.h"
#include "EnergyArray.h"
//#incl ude "CpptrajStdio.h" // DEBUG

void PotentialTerm_Bond::addBonds(BondArray const& bonds, CharMask const& maskIn) {
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    if (maskIn.AtomInCharMask( bnd->A1() ) ||
        maskIn.AtomInCharMask( bnd->A2() ))
    {
      //mprintf("DEBUG: Bond %i to %i\n", bnd->A1()+1, bnd->A2()+1);
      activeBonds_.push_back( *bnd );
    }
  }
}

/** Set up Hooke's law bond term. */
int PotentialTerm_Bond::SetupTerm(Topology const& topIn, Box const& boxIn,
                                  CharMask const& maskIn, EnergyArray& Earray)
{
  activeBonds_.clear();
  addBonds( topIn.Bonds(),  maskIn );
  addBonds( topIn.BondsH(), maskIn );

  bondParm_ = &(topIn.BondParm());
  Ebond_ = Earray.AddType( EnergyArray::E_BOND );

  return 0;
}

/** Calculate Hooke's law bond force. */
void PotentialTerm_Bond::CalcForce(Frame& frameIn, CharMask const& maskIn) const {
  *Ebond_ = 0.0;
  for (BondArray::const_iterator bnd = activeBonds_.begin(); bnd != activeBonds_.end(); ++bnd)
  {
    BondParmType BP = (*bondParm_)[ bnd->Idx() ];
    //Vec3 const& XYZ0 = Xarray[ bnd->A1() ];
    //Vec3 const& XYZ1 = Xarray[ bnd->A2() ];
    const double* XYZ0 = frameIn.XYZ( bnd->A1() );
    const double* XYZ1 = frameIn.XYZ( bnd->A2() );
    double rx = XYZ0[0] - XYZ1[0];
    double ry = XYZ0[1] - XYZ1[1];
    double rz = XYZ0[2] - XYZ1[2];
    double r2 = rx*rx + ry*ry + rz*rz;
    if (r2 > 0.0) {
      double r2inv = 1.0/r2;
      double r = sqrt(r2);
      //mprintf("DBG: %i A1=%i A2=%i R=%g\n", iteration, bnd->A1()+1, bnd->A2()+1, r);
      double rinv = r * r2inv;

      double db = r - BP.Req();
      double df = BP.Rk() * db;
      double e = df * db;
      *Ebond_ += e;

      df *= 2.0 * rinv;

      double dfx = df * rx;
      double dfy = df * ry;
      double dfz = df * rz;

      if (maskIn.AtomInCharMask(bnd->A1())) {
        double* fxyz = frameIn.fAddress() + (3*bnd->A1());
        fxyz[0] -= dfx;
        fxyz[1] -= dfy;
        fxyz[2] -= dfz;
      }
     if (maskIn.AtomInCharMask(bnd->A2())) {
        double* fxyz = frameIn.fAddress() + (3*bnd->A2());
        fxyz[0] += dfx;
        fxyz[1] += dfy;
        fxyz[2] += dfz;
      }
    }
  }
}
