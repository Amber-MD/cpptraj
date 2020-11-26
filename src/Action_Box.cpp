#include "Action_Box.h"
#include "CpptrajStdio.h"

Action_Box::Action_Box() :
  mode_(SET)
{
  for (int i = 0; i < 6; i++) {
    xyzabg_[i] = 0;
    setVar_[i] = false;
  }
}

void Action_Box::Help() const {
  mprintf("\t{[x <xval>] [y <yval>] [z <zval>] {[alpha <a>] [beta <b>] [gamma <g>]\n"
          "\t [truncoct]} | nobox | auto [offset <offset>] [radii {vdw|gb|parse|none}]}\n"
          "  For each input frame, replace any box information with the information given.\n"
          "  If 'truncoct' is specified, alpha, beta, and gamma will be set to the\n"
          "  appropriate angle for a truncated octahedral box. If 'nobox' is specified,\n"
          "  all existing box information will be removed. If 'auto' is specified, an\n"
          "  orthogonal box will be set for existing atoms using the specified distance\n"
          "  offset value, ensuring specified radii (default vdw) are enclosed.\n");
}

// Action_Box::Init()
Action::RetType Action_Box::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  if ( actionArgs.hasKey("nobox") )
    mode_ = REMOVE;
  else if (actionArgs.hasKey("auto")) {
    offset_ = actionArgs.getKeyDouble("offset", 0.0);
    if (offset_ < 0) {
      mprinterr("Error: Offset for auto must be >= 0.\n");
      return Action::ERR;
    }
    mode_ = AUTO;
    radiiMode_ = UNSPECIFIED;
    // NOTE: Set angles to 90 here so we can set parm box type to
    //       ORTHO in Setup().
    xyzabg_[3] = 90.0;
    xyzabg_[4] = 90.0;
    xyzabg_[5] = 90.0;
    xyzabg_[0] = 1.0;
    xyzabg_[1] = 1.0;
    xyzabg_[2] = 1.0;
    std::string rstr = actionArgs.GetStringKey("radii");
    if (!rstr.empty()) {
      if (rstr == "vdw")
        radiiMode_ = VDW;
      else if (rstr == "parse")
        radiiMode_ = PARSE;
      else if (rstr == "gb")
        radiiMode_ = GB;
      else if (rstr == "none")
        radiiMode_ = NONE;
      else {
        mprinterr("Error: Unrecognized radii type: %s\n", rstr.c_str());
        return Action::ERR;
      }
    }
  } else {
    mode_ = SET;
    // TODO check for bad args?
    if (actionArgs.Contains("x")) { xyzabg_[0] = actionArgs.getKeyDouble("x", 0.0); setVar_[0] = true; }
    if (actionArgs.Contains("y")) { xyzabg_[1] = actionArgs.getKeyDouble("y", 0.0); setVar_[1] = true; }
    if (actionArgs.Contains("z")) { xyzabg_[2] = actionArgs.getKeyDouble("z", 0.0); setVar_[2] = true; }
    if (actionArgs.Contains("alpha")) { xyzabg_[3] = actionArgs.getKeyDouble("alpha", 0.0); setVar_[3] = true; }
    if (actionArgs.Contains("beta"))  { xyzabg_[4] = actionArgs.getKeyDouble("beta",  0.0); setVar_[4] = true; }
    if (actionArgs.Contains("gamma")) { xyzabg_[5] = actionArgs.getKeyDouble("gamma", 0.0); setVar_[5] = true; }
    if (actionArgs.hasKey("truncoct")) {
      xyzabg_[3] = Box::TruncatedOctAngle();
      xyzabg_[4] = xyzabg_[3];
      xyzabg_[5] = xyzabg_[3];
      setVar_[3] = true;
      setVar_[4] = true;
      setVar_[5] = true;
      // All lengths need to be the same
      if (setVar_[1]) mprintf("Warning: Only 'x' used for 'truncoct'\n");
      if (setVar_[2]) mprintf("Warning: Only 'x' used for 'truncoct'\n");
      if (setVar_[0]) {
        xyzabg_[1] = xyzabg_[0];
        xyzabg_[2] = xyzabg_[0];
      }
      setVar_[1] = false;
      setVar_[2] = false;
    }
  }

  mprintf("    BOX:");
  if (mode_ == REMOVE)
    mprintf(" Removing box information.\n");
  else if (mode_ == AUTO) {
    mprintf(" Setting orthogonal box for atoms using offset of %g Ang\n", offset_);
    switch (radiiMode_) {
      case GB    : mprintf("\tUsing GB radii.\n"); break;
      case PARSE : mprintf("\tUsing PARSE radii.\n"); break;
      case VDW   : mprintf("\tUsing VDW radii.\n"); break;
      case NONE  : mprintf("\tNot using atomic radii.\n"); break;
      case UNSPECIFIED:
        mprintf("\tWill use VDW, GB, or PARSE radii if available (with that priority).\n");
        break;
    }
  } else {
    if (setVar_[0]) mprintf(" X=%.3f", xyzabg_[0]);
    if (setVar_[1]) mprintf(" Y=%.3f", xyzabg_[1]);
    if (setVar_[2]) mprintf(" Z=%.3f", xyzabg_[2]);
    if (setVar_[3]) mprintf(" A=%.3f", xyzabg_[3]);
    if (setVar_[4]) mprintf(" B=%.3f", xyzabg_[4]);
    if (setVar_[5]) mprintf(" G=%.3f", xyzabg_[5]);
    mprintf("\n");
  }
  return Action::OK;
}

/** Set missing box information from incoming box. */
void Action_Box::SetMissingInfo(Box const& boxIn) {
  if (!setVar_[0]) xyzabg_[0] = boxIn.BoxX();
  if (!setVar_[1]) xyzabg_[1] = boxIn.BoxY();
  if (!setVar_[2]) xyzabg_[2] = boxIn.BoxZ();
  if (!setVar_[3]) xyzabg_[3] = boxIn.Alpha();
  if (!setVar_[4]) xyzabg_[4] = boxIn.Beta();
  if (!setVar_[5]) xyzabg_[5] = boxIn.Gamma();
}

// Action_Box::Setup()
Action::RetType Action_Box::Setup(ActionSetup& setup) {
  cInfo_ = setup.CoordInfo();
  if (mode_ == REMOVE) {
    mprintf("\tRemoving box info.\n");
    cInfo_.SetBox( Box() );
  } else {
    // SET, AUTO
    // Fill in missing box information from current box
    SetMissingInfo( setup.CoordInfo().TrajBox() );
    Box pbox;
    pbox.SetupFromXyzAbg( xyzabg_ );
    mprintf("\tNew box type is %s\n", pbox.TypeName() );
    cInfo_.SetBox( pbox );
    // Get radii for AUTO
    if (mode_ == AUTO) {
      if (radiiMode_ == VDW && !setup.Top().Nonbond().HasNonbond()) {
        mprintf("Warning: No VDW radii in topology %s; skipping.\n", setup.Top().c_str());
        return Action::SKIP;
      }
      RadiiType modeToUse = radiiMode_;
      if (modeToUse == UNSPECIFIED) {
        // If VDW radii present, use those.
        if (setup.Top().Nonbond().HasNonbond())
          modeToUse = VDW;
        else if (setup.Top().Natom() > 0 && setup.Top()[0].GBRadius() > 0)
          modeToUse = GB;
        else
          modeToUse = PARSE;
      }
      switch (modeToUse) {
        case GB    : mprintf("\tUsing GB radii.\n"); break;
        case PARSE : mprintf("\tUsing PARSE radii.\n"); break;
        case VDW   : mprintf("\tUsing VDW radii.\n"); break;
        case UNSPECIFIED:
        case NONE:
          mprintf("\tNot using atomic radii.\n");
          break;
      }
      Radii_.clear();
      Radii_.reserve( setup.Top().Natom() );
      for (int atnum = 0; atnum != setup.Top().Natom(); ++atnum) {
        switch (modeToUse) {
          case GB   : Radii_.push_back( setup.Top()[atnum].GBRadius()    ); break;
          case PARSE: Radii_.push_back( setup.Top()[atnum].ParseRadius() ); break;
          case VDW  : Radii_.push_back( setup.Top().GetVDWradius(atnum)  ); break;
          case UNSPECIFIED:
          case NONE:
            Radii_.push_back( 0.0 );
            break;
        }
      }
    } 
  }
  setup.SetCoordInfo( &cInfo_ );
  return Action::MODIFY_TOPOLOGY;
}

// Action_Box::DoAction()
Action::RetType Action_Box::DoAction(int frameNum, ActionFrame& frm) {
  if (mode_ == REMOVE) {
    frm.ModifyFrm().SetBox( Box() );
  } else if (mode_ == AUTO) {
    frm.ModifyFrm().SetOrthoBoundingBox(Radii_, offset_);
  } else {
    SetMissingInfo( frm.Frm().BoxCrd() );
    frm.ModifyFrm().ModifyBox().SetupFromXyzAbg( xyzabg_ );
  }
  return Action::MODIFY_COORDS;
}
