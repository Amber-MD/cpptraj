#include <algorithm> // sort
#include <cmath> // fabs
#include "Constraints.h"
#include "CpptrajStdio.h"
#include "Constants.h"

/// CONSTRUCTOR
Constraints::Constraints() :
  dt_(0.002),
  epsilon_(0.0000001),
  EPS_(epsilon_ / (dt_ * Constants::AMBERTIME_TO_PS)),
  shakeType_(OFF),
  degrees_of_freedom_(0)
{}

static const char* shakeString_[] = {
  "off", "bonds to H", "all bonds"
};

const char* Constraints::shakeString() const { return shakeString_[shakeType_]; }

const char* Constraints::constraintArgs = "[ntc <#>]";

// Constraints::InitConstraints()
int Constraints::InitConstraints(ArgList& argIn) {
  int ntc = argIn.getKeyInt("ntc",-1);
  if (ntc != -1) {
    if (ntc < 1 || ntc > 3) {
      mprinterr("Error: ntc must be 1 (off), 2 (bonds to H), or 3 (all bonds).\n");
      return 1;
    }
    shakeType_ = (ShakeType)(ntc - 1);
  } else
    shakeType_ = OFF;

  return 0;
}

const char* Constraints::rattleArgs = "[dt <time>] [epsilon <eps>]";

//  Constraints::InitRattle()
int Constraints::InitRattle(ArgList& argIn) {
  dt_ = argIn.getKeyDouble("dt", 0.002);
  // FIXME - check this conversion
  double dtx = dt_ * Constants::AMBERTIME_TO_PS;
  epsilon_ = argIn.getKeyDouble("epsilon", 0.0000001);
  EPS_ = epsilon_ / dtx;
  return 0;
}

// Constraints::AddBonds()
int Constraints::AddBonds(BondArray const& bonds, Topology const& Top,
                          CharMask const& cMask)
{
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    if (cMask.AtomInCharMask(bnd->A1()) && cMask.AtomInCharMask(bnd->A2()))
    {
      int idx = bnd->Idx();
      if (idx < 0) {
        mprinterr("Error: No bond parameters for atoms %s and %s\n",
                  Top.TruncAtomNameNum(bnd->A1()).c_str(),
                  Top.TruncAtomNameNum(bnd->A2()).c_str());
        return 1;
      }
      double req = Top.BondParm()[idx].Req();
      if (req < Constants::SMALL) {
        mprinterr("Error: Zero bond length for %s to %s\n",
                  Top.TruncAtomNameNum(bnd->A1()).c_str(),
                  Top.TruncAtomNameNum(bnd->A2()).c_str());
        return 1;
      }
      Bonds_.push_back( Cbond(bnd->A1(), bnd->A2(), req) );
    }
  }
  return 0;
}

// Constraints::SetupConstraints()
int Constraints::SetupConstraints(AtomMask const& mask, Topology const& top)
{
  Bonds_.clear();
  unsigned int constrained_bonds_to_h = 0;
  unsigned int constrained_heavy_bonds = 0;
  if (shakeType_ > OFF) {
    // Both atoms must be in the mask. TODO should be either? TODO Warn about partial bonds?
    CharMask cMask(mask.ConvertToCharMask(), mask.Nselected());
    if (AddBonds(top.BondsH(), top, cMask)) return 1;
    constrained_bonds_to_h = Bonds_.size();
    mprintf("\tConstraints on %u bonds to hydrogen", constrained_bonds_to_h);
    if (shakeType_ > BONDS_TO_H) {
      if (AddBonds(top.Bonds(), top, cMask)) return 1;
      constrained_heavy_bonds = Bonds_.size() - constrained_bonds_to_h;
      std::sort( Bonds_.begin(), Bonds_.end() );
      mprintf(", %u heavy atom bonds, %zu bonds total", constrained_heavy_bonds, Bonds_.size());
    }
    mprintf(".\n");
  }
  // Calculate degrees of freedom
  degrees_of_freedom_ = (3 * mask.Nselected()) - 
                        constrained_bonds_to_h -
                        constrained_heavy_bonds;
    mprintf("\t# of degrees of freedom = %i\n", degrees_of_freedom_);
  return 0;
}

// Constraints::Rattle2()
int Constraints::Rattle2(Frame& frameIn)
const
{
  static const double FAC = 1.2;
  // Main loop
  bool done = false;
  int niterations = 0;
  static const int maxIterations = 1000;
  double* Vel = frameIn.vAddress();
  while (!done && (niterations < maxIterations))
  {
    niterations++;
    done = true;
    // Loop over selected constrained bonds
    for (Carray::const_iterator bnd = Bonds_.begin(); bnd != Bonds_.end(); ++bnd)
    {
        // Get the bond equilibrium length
        double req = bnd->req_;
        // Calculate bond vector
        Vec3 dR = Vec3(frameIn.XYZ( bnd->at2_)) - Vec3(frameIn.XYZ( bnd->at1_ ));
        // Get velocities for atom 1 (A)
        double* VA = Vel + (bnd->at1_*3);
        // Get velocities for atom 2 (B)
        double* VB = Vel + (bnd->at2_*3);
        // Calculate delta V
        Vec3 dV = Vec3(VB) - Vec3(VA);
        double dot = dR * dV;
        double rma = 1.0 / frameIn.Mass( bnd->at1_ );
        double rmb = 1.0 / frameIn.Mass( bnd->at2_ );
        double term = -dot * FAC / ((rma + rmb) * req * req);
        if (fabs(term) > EPS_) {
          done = false;
          Vec3 dTerm = dR * term;

          VA[0] -= dTerm[0] * rma;
          VA[1] -= dTerm[1] * rma;
          VA[2] -= dTerm[2] * rma;

          VB[0] += dTerm[0] * rmb; 
          VB[1] += dTerm[1] * rmb; 
          VB[2] += dTerm[2] * rmb;
        }
    } // END loop over bonds
  } // END main rattle loop
  if (niterations > maxIterations) {
    mprinterr("Error: RATTLE took more than %i iterations.\n", maxIterations);
    return 1;
  }
  return 0;
}
