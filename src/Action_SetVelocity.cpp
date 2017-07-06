#include <cmath> // sqrt
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
  double dt = actionArgs.getKeyDouble("dt", 0.002);
  // FIXME
  dt *= Constants::AMBERTIME_TO_PS; // dtx
  double epsilon = actionArgs.getKeyDouble("epsilon", 0.0000001);
  EPS_ = epsilon / dt;
  // Masks
  Mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    SETVELOCITY: Assigning velocities for atoms in mask '%s'\n", Mask_.MaskString());
  mprintf("\tTemperature= %.2f, using Maxwellian distribution.\n", tempi_);
  if (ig_ != -1)
    mprintf("\tRandom seed is %i\n", ig_);
  return Action::OK;
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
  cMask_ = CharMask(Mask_.ConvertToCharMask(), Mask_.Nselected());
  SD_.clear();
  SD_.reserve( Mask_.Nselected() );
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
  // Always add velocity info even if not strictly necessary
  cInfo_ = setup.CoordInfo();
  cInfo_.SetVelocity( true );
  newFrame_.SetupFrameV( setup.Top().Atoms(), cInfo_ );
  setup.SetCoordInfo( &cInfo_ );
  return Action::MODIFY_TOPOLOGY;
}

int Action_SetVelocity::Rattle2(Frame& frameIn, BondParmArray const& parmIn, BondArray const& bonds)
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
    // Loop over bonds
    for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
    {
      // Both atoms must be selected
      if (cMask_.AtomInCharMask(bnd->A1()) && cMask_.AtomInCharMask(bnd->A2()))
      {
        // Get the force constant. Assume they were checked for in Setup()
        double rk = parmIn[bnd->Idx()].Rk();
        // Calculate bond vector
        Vec3 dR = Vec3(frameIn.XYZ( bnd->A2())) - Vec3(frameIn.XYZ( bnd->A1()));
        // Get velocities for atom 1 (A)
        double* VA = Vel + (bnd->A1()*3);
        // Get velocities for atom 2 (B)
        double* VB = Vel + (bnd->A2()*3);
        // Calculate delta V
        Vec3 dV = Vec3(VB) - Vec3(VA);
        double dot = dR * dV;
        double rma = 1.0 / frameIn.Mass(bnd->A1());
        double rmb = 1.0 / frameIn.Mass(bnd->A2());
        double term = -dot * FAC / ((rma + rmb) * rk * rk);
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
      } // END both atoms selected 
    } // END loop over bonds
  } // END main rattle loop
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
    std::vector<double>::const_iterator sd = SD_.begin(); 
    for (AtomMask::const_iterator atom = Mask_.begin(); atom != Mask_.end(); ++atom, ++sd)
    {
      double* V = newFrame_.vAddress() + (*atom * 3);
      V[0] = RN_.rn_gauss(0.0, *sd);
      V[1] = RN_.rn_gauss(0.0, *sd);
      V[2] = RN_.rn_gauss(0.0, *sd);
    }
  }
  frm.SetFrame( &newFrame_ );
  return Action::MODIFY_COORDS;
}
