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
          "\t getbox {ucell|frac|shape} [name <setname>] [out <file>]\n"
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
  std::string dsname;
  DataFile* outfile = 0;
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
    mode_ = GETBOX;
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
    dsname = actionArgs.GetStringKey("name");
    outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  } else {
    mode_ = SET;
    // TODO check for bad args?
    if (boxArgs_.SetBoxArgs( actionArgs )) return Action::ERR;
  }

  // Create set
  if (mode_ == GETBOX) {
    set_ = init.DSL().AddSet( DataSet::MAT3X3, dsname, "BOX" );
    if (set_ == 0) {
      mprinterr("Error: Could not create box set.\n");
      return Action::ERR;
    }
    if (outfile != 0) outfile->AddDataSet( set_ );
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
  } else if (mode_ == GETBOX) {
    static const char* getstr[] = { "unit cell", "fractional cell", "shape matrix" };
    mprintf(" Getting %s information from box.\n", getstr[getmode_]);
    mprintf("\tOutput set: %s\n", set_->legend());
    if (outfile != 0) mprintf("\tOutput file: %s\n", outfile->DataFilename().full());
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
  } else if (mode_ == GETBOX) {
  // Check box type
    if (!setup.CoordInfo().TrajBox().HasBox()) {
      mprintf("Warning: Topology %s does not contain box information.\n",
              setup.Top().c_str());
      return Action::SKIP;
    }
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
  } else if (mode_ == GETBOX) {
    if (getmode_ == GET_UNITCELL) {
      Matrix_3x3 const& ucell = frm.Frm().BoxCrd().UnitCell();
      set_->Add( frameNum, ucell.Dptr() );
    } else if (getmode_ == GET_FRACCELL) {
      Matrix_3x3 const& frac = frm.Frm().BoxCrd().FracCell();
      set_->Add( frameNum, frac.Dptr() );
    } else if (getmode_ == GET_SHAPE) {
      double shape[6];
      frm.Frm().BoxCrd().GetSymmetricShapeMatrix(shape);
      double ucell[9];
      ucell[0] = shape[0]; ucell[1] = shape[1]; ucell[2] = shape[3];
      ucell[3] = shape[1]; ucell[4] = shape[2]; ucell[5] = shape[4];
      ucell[6] = shape[3]; ucell[7] = shape[4]; ucell[8] = shape[5];
      set_->Add( frameNum, ucell );
    }
  } else {
    boxArgs_.SetMissingInfo( frm.Frm().BoxCrd() );
    frm.ModifyFrm().ModifyBox().AssignFromXyzAbg( boxArgs_.XyzAbg() );
  }
  return Action::MODIFY_COORDS;
}
