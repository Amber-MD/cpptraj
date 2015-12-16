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
