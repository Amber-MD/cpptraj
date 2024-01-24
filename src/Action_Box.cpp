#include "Action_Box.h"
#include "CpptrajStdio.h"

Action_Box::Action_Box() :
  mode_(SET),
  set_(0),
  getmode_(GET_UNITCELL)
{}

void Action_Box::Help() const {
  mprintf("\t{%s |\n"
          "\t %s |\n"
          "\t nobox |\n"
          "\t auto [offset <offset>] [radii {vdw|gb|parse|none}] |\n"
          "\t getbox {ucell|frac|shape}\n"
          "\t}\n",
          BoxArgs::Keywords_XyzAbg(), BoxArgs::Keywords_TruncOct());
  mprintf("  For each input frame, replace any box information with the information given.\n"
          "  If 'truncoct' is specified, alpha, beta, and gamma will be set to the\n"
          "  appropriate angle for a truncated octahedral box.\n"
          "  If 'nobox' is specified, all existing box information will be removed.\n"
          "  If 'auto' is specified, an orthogonal box will be set for existing atoms\n"
          "  using the specified distance offset value, ensuring specified radii\n"
          "  (default vdw) are enclosed.\n"
          "  If 'getbox' is specified, the existing box information will be saved to\n"
          "  a data set.\n");
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
    boxArgs_.SetAngles( 90.0 );
    boxArgs_.SetLengths( 1.0 );
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
  } else if (actionArgs.Contains("getbox")) {
    std::string getbox = actionArgs.GetStringKey("getbox");
    if (getbox == "ucell")
      getmode_ = GET_UNITCELL;
    else if (getbox == "frac")
      getmode_ = GET_FRACCELL;
    else if (getbox == "shape")
      getmode_ = GET_SHAPE;
    else {
      mprinterr("Error: Expected getbox {ucell|frac|shape}, got %s\n", getbox.c_str());
      return Action::ERR;
    }
  } else {
    mode_ = SET;
    // TODO check for bad args?
    if (boxArgs_.SetBoxArgs( actionArgs )) return Action::ERR;
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
    boxArgs_.PrintXyzAbg();
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
    // Fill in missing box information from current box
    if (boxArgs_.SetMissingInfo( setup.CoordInfo().TrajBox() ))
      return Action::ERR;
    Box pbox;
    pbox.SetupFromXyzAbg( boxArgs_.XyzAbg() );
    mprintf("\tNew box type is %s\n", pbox.CellShapeName() );
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
    frm.ModifyFrm().ModifyBox().SetNoBox();
  } else if (mode_ == AUTO) {
    frm.ModifyFrm().SetOrthoBoundingBox(Radii_, offset_);
  } else {
    boxArgs_.SetMissingInfo( frm.Frm().BoxCrd() );
    frm.ModifyFrm().ModifyBox().AssignFromXyzAbg( boxArgs_.XyzAbg() );
  }
  return Action::MODIFY_COORDS;
}
