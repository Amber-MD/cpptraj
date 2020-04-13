#include "Action_Box.h"
#include "CpptrajStdio.h"

Action_Box::Action_Box() : mode_(SET) {}

void Action_Box::Help() const {
  mprintf("\t{[x <xval>] [y <yval>] [z <zval>] {[alpha <a>] [beta <b>] [gamma <g>]\n"
          "\t [truncoct]} | nobox | auto <offset>}\n"
          "  For each input frame, replace any box information with the information given.\n"
          "  If 'truncoct' is specified, alpha, beta, and gamma will be set to the\n"
          "  appropriate angle for a truncated octahedral box. If 'nobox' is specified,\n"
          "  all existing box information will be removed. If 'auto' is specified, an\n"
          "  orthogonal box will be set for existing atoms using the specified distance\n"
          "  offset value.\n");
}

// Action_Box::Init()
Action::RetType Action_Box::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  if ( actionArgs.hasKey("nobox") )
    mode_ = REMOVE;
  else if (actionArgs.Contains("auto")) {
    offset_ = actionArgs.getKeyDouble("auto", -1.0);
    if (offset_ < 0) {
      mprinterr("Error: Offset for auto must be >= 0.\n");
      return Action::ERR;
    }
    mode_ = AUTO;
    box_.SetAlpha(90.0);
    box_.SetBeta(90.0);
    box_.SetGamma(90.0);
    box_.SetX(1.0);
    box_.SetY(1.0);
    box_.SetZ(1.0);
  } else {
    mode_ = SET;
    box_.SetX( actionArgs.getKeyDouble("x", 0.0) );
    box_.SetY( actionArgs.getKeyDouble("y", 0.0) );
    box_.SetZ( actionArgs.getKeyDouble("z", 0.0) );
    box_.SetAlpha( actionArgs.getKeyDouble("alpha", 0.0) );
    box_.SetBeta(  actionArgs.getKeyDouble("beta",  0.0) );
    box_.SetGamma( actionArgs.getKeyDouble("gamma", 0.0) );
    if (actionArgs.hasKey("truncoct")) box_.SetTruncOct();
  }

  mprintf("    BOX:");
  if (mode_ == REMOVE)
    mprintf(" Removing box information.\n");
  else if (mode_ == AUTO)
    mprintf(" Setting orthogonal box for atoms using offset of %g Ang\n", offset_);
  else {
    if (box_.BoxX() > 0) mprintf(" X=%.3f", box_.BoxX());
    if (box_.BoxY() > 0) mprintf(" Y=%.3f", box_.BoxY());
    if (box_.BoxZ() > 0) mprintf(" Z=%.3f", box_.BoxZ());
    if (box_.Alpha() > 0) mprintf(" A=%.3f", box_.Alpha());
    if (box_.Beta() > 0) mprintf(" B=%.3f", box_.Beta());
    if (box_.Gamma() > 0) mprintf(" G=%.3f", box_.Gamma());
    mprintf("\n");
  }
  return Action::OK;
}

// Action_Box::Setup()
Action::RetType Action_Box::Setup(ActionSetup& setup) {
  cInfo_ = setup.CoordInfo();
  if (mode_ == REMOVE) {
    mprintf("\tRemoving box info.\n");
    cInfo_.SetBox( Box() );
  } else {
    // SET, AUTO
    Box pbox( box_ );
    // Fill in missing box information from current box 
    pbox.SetMissingInfo( setup.CoordInfo().TrajBox() );
    mprintf("\tNew box type is %s\n", pbox.TypeName() );
    cInfo_.SetBox( pbox );
  }
  setup.SetCoordInfo( &cInfo_ );
  return Action::MODIFY_TOPOLOGY;
}

// Action_Box::DoAction()
Action::RetType Action_Box::DoAction(int frameNum, ActionFrame& frm) {
  if (mode_ == REMOVE) {
    frm.ModifyFrm().SetBox( Box() );
  } else if (mode_ == AUTO) {
    Box fbox( box_ );
    int atom = 0;
    Vec3 min(frm.Frm().XYZ( atom ));
    Vec3 max(min);
    for (; atom != frm.Frm().Natom(); ++atom)
    {
      const double* xyz = frm.Frm().XYZ( atom );
      if (xyz[0] < min[0]) min[0] = xyz[0];
      if (xyz[0] > max[0]) max[0] = xyz[0];
      if (xyz[1] < min[1]) min[1] = xyz[1];
      if (xyz[1] > max[1]) max[1] = xyz[1];
      if (xyz[2] < min[2]) min[2] = xyz[2];
      if (xyz[2] > max[2]) max[2] = xyz[2];
    }
    min -= offset_;
    max += offset_;
    fbox.SetX(max[0] - min[0]);
    fbox.SetY(max[1] - min[1]);
    fbox.SetZ(max[2] - min[2]);
    frm.ModifyFrm().SetBox( fbox );
  } else {
    Box fbox( box_ );
    fbox.SetMissingInfo( frm.Frm().BoxCrd() );
    frm.ModifyFrm().SetBox( fbox );
  }
  return Action::MODIFY_COORDS;
}
