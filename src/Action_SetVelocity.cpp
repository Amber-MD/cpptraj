#include <cmath> // sqrt
#include <algorithm> // sort
#include "Action_SetVelocity.h"
#include "CpptrajStdio.h"
#include "Constants.h"

Action_SetVelocity::Action_SetVelocity() : tempi_(0.0) {}

void Action_SetVelocity::Help() const {
  mprintf("\t[<mask>] [tempi <temperature>] [ig <random seed>]\n"
          "  Set velocities in frame for atoms in <mask> using Maxwellian distribution\n" 
          "  based on given temperature.\n");
}

// Action_SetVelocity::Init()
Action::RetType Action_SetVelocity::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Keywords
  tempi_ = actionArgs.getKeyDouble("tempi", 300.0);
  int ig_ = actionArgs.getKeyInt("ig", -1);
  RN_.rn_set( ig_ );
  ntc_ = actionArgs.getKeyInt("ntc", 1);
  if (ntc_ < 1 || ntc_ > 3) {
    mprinterr("Error: 'ntc' must be between 1 and 3.\n");
    return Action::ERR;
  }
  // If bonds will be constrained, get some more info.
  if (ntc_ > 1) {
    double dt = actionArgs.getKeyDouble("dt", 0.002);
    // FIXME - check this conversion
    dt *= Constants::AMBERTIME_TO_PS; // dtx
    double epsilon = actionArgs.getKeyDouble("epsilon", 0.0000001);
    EPS_ = epsilon / dt;
  } else
    EPS_ = 1.0;
  // Masks
  Mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    SETVELOCITY: Assigning velocities for atoms in mask '%s'\n", Mask_.MaskString());
  mprintf("\tTemperature= %.2f, using Maxwellian distribution.\n", tempi_);
  if (ig_ != -1)
    mprintf("\tRandom seed is %i\n", ig_);
  return Action::OK;
}

//  Action_SetVelocity::AddBonds()
int Action_SetVelocity::AddBonds(BondArray const& bonds, Topology const& Top,
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
        return Action::ERR;
      }
      Bonds_.push_back( Cbond(bnd->A1(), bnd->A2(), req) );
    }
  }
  return 0;
}

// Action_SetVelocity::Setup()
Action::RetType Action_SetVelocity::Setup(ActionSetup& setup) {
  // Masks
  if (setup.Top().SetupIntegerMask( Mask_ )) return Action::ERR;
  Mask_.MaskInfo();
  if (Mask_.None()) {
    mprintf("Warning: No atoms selected in [%s]\n", Mask_.MaskString());
    return Action::SKIP;
  }
  SD_.clear();
  SD_.reserve( Mask_.Nselected() );
  // Store mass-related values for atoms
  double boltz = Constants::GASK_KCAL * tempi_;
  for (AtomMask::const_iterator atom = Mask_.begin(); atom != Mask_.end(); ++atom)
  {
    double mass_inv;
    double mass = setup.Top()[*atom].Mass();
    if ( mass < Constants::SMALL )
      mass_inv = 0.0;
    else
      mass_inv = 1.0 / mass;
    SD_.push_back( sqrt(boltz * mass_inv) );
  }
  // Save bond info if using constraints
  if (ntc_ > 1) {
    Bonds_.clear();
    // Both bonds must be in the mask. TODO warn about partial bonds?
    CharMask cMask(Mask_.ConvertToCharMask(), Mask_.Nselected());
    if (AddBonds(setup.Top().BondsH(), setup.Top(), cMask)) return Action::ERR;
    mprintf("\tConstraints on %zu bonds to hydrogen", Bonds_.size());
    if (ntc_ > 2) {
      if (AddBonds(setup.Top().Bonds(), setup.Top(), cMask)) return Action::ERR;
      std::sort( Bonds_.begin(), Bonds_.end() );
      mprintf(", %zu bonds total", Bonds_.size());
    }
    mprintf(".\n");
  }
  // Always add velocity info even if not strictly necessary
  cInfo_ = setup.CoordInfo();
  cInfo_.SetVelocity( true );
  newFrame_.SetupFrameV( setup.Top().Atoms(), cInfo_ );
  setup.SetCoordInfo( &cInfo_ );
  return Action::MODIFY_TOPOLOGY;
}

// Action_SetVelocity::Rattle2()
int Action_SetVelocity::Rattle2(Frame& frameIn)
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

// Action_SetVelocity::DoAction()
Action::RetType Action_SetVelocity::DoAction(int frameNum, ActionFrame& frm) {
  std::copy( frm.Frm().xAddress(), frm.Frm().xAddress() + frm.Frm().size(), newFrame_.xAddress() );
  if (frm.Frm().HasVelocity())
    std::copy( frm.Frm().vAddress(), frm.Frm().vAddress() + frm.Frm().size(),
               newFrame_.vAddress() );
  if (tempi_ < Constants::SMALL) {
    for (AtomMask::const_iterator atom = Mask_.begin(); atom != Mask_.end(); ++atom)
    {
      double* V = newFrame_.vAddress() + (*atom * 3);
      V[0] = 0.0;
      V[1] = 0.0;
      V[2] = 0.0;
    }
  } else {
    Darray::const_iterator sd = SD_.begin(); 
    for (AtomMask::const_iterator atom = Mask_.begin(); atom != Mask_.end(); ++atom, ++sd)
    {
      double* V = newFrame_.vAddress() + (*atom * 3);
      V[0] = RN_.rn_gauss(0.0, *sd);
      V[1] = RN_.rn_gauss(0.0, *sd);
      V[2] = RN_.rn_gauss(0.0, *sd);
    }
  }
  if (ntc_ > 1)
    Rattle2( newFrame_ );
  frm.SetFrame( &newFrame_ );
  return Action::MODIFY_COORDS;
}
