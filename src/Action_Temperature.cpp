#include "Action_Temperature.h"
#include "CpptrajStdio.h"

Action_Temperature::Action_Temperature() :
  Tdata_(0),
  mode_(CALC_ONLY),
  dof_offset_(0),
  removeTrans_(false),
  removeRot_(false)
{}

void Action_Temperature::Help() const {
  mprintf("\t[<name>] [out <filename>]\n"
          "\t{ frame |\n"
          "\t  [<mask>] %s [update] [remove {trans|rot|both}]\n"
          "\t}\n", Constraints::constraintArgs);
  mprintf("  Calculate temperature in frame based on velocity information. If\n"
          "  'update' is specified, update frame temperature too. If 'frame'\n"
          "  is specified just use frame temperature (e.g. read in from a\n"
          "  REMD trajectory).\n"
          "  The 'ntc' keyword can be used to correct\n"
          "  for lost degrees of freedom due to SHAKE constraints (2 = bonds to\n"
          "  hydrogen, 3 = all bonds). The 'remove' keyword can be used to\n"
          "  account for removed translational and/or rotational degrees of\n"
          "  freedom.\n");
}

// Action_Temperature::Init()
Action::RetType Action_Temperature::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Keywords
  if (actionArgs.hasKey("frame"))
    mode_ = FROM_FRAME;
  else {
    mode_ = CALC_ONLY;
    if (cons_.InitConstraints( actionArgs )) return Action::ERR;
  }
  if (mode_ == CALC_ONLY && actionArgs.hasKey("update"))
    mode_ = CALC_AND_MODIFY;
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );

  // Take into account removal of global degrees of freedom?
  removeTrans_ = false;
  removeRot_ = false;
  if (mode_ != FROM_FRAME) {
    std::string removearg = actionArgs.GetStringKey("remove");
    if (!removearg.empty()) {
      if (removearg == "trans")
        removeTrans_ = true;
      else if (removearg == "rot")
        removeRot_ = true;
      else if (removearg == "both") {
        removeTrans_ = true;
        removeRot_ = true;
      } else {
        mprinterr("Error: Unrecognized arg for 'remove' keyword: %s\n", removearg.c_str());
        return Action::ERR;
      }
    }
    // Masks
    if (Mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  }
  // DataSet 
  Tdata_ =  init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "Tdata");
  if (Tdata_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddDataSet( Tdata_ );
  
  if (mode_ == FROM_FRAME) {
    mprintf("    TEMPERATURE: Frame temperatures will be saved in data set %s\n",
             Tdata_->legend());
  } else {
    mprintf("    TEMPERATURE: Calculate temperature for atoms in mask [%s]\n", Mask_.MaskString());
    if (mode_ == CALC_AND_MODIFY)
      mprintf("\tAny existing temperature in Frames will be overwritten.\n");
    mprintf("\tConstraints: %s\n", cons_.shakeString());
    if (removeTrans_) mprintf("\tAssuming translational degs. of freedom removed.\n");
    if (removeRot_) mprintf("\tAssuming rotational degs. of freedom removed.\n");
  }
  return Action::OK;
}

// Action_Temperature::Setup()
Action::RetType Action_Temperature::Setup(ActionSetup& setup) {
  cInfo_ = setup.CoordInfo();
  Action::RetType ret = Action::OK;
  if (mode_ == FROM_FRAME) {
    // Get temperature from frame.
    if (!cInfo_.HasTemp()) {
      mprintf("Warning: No temperature information in Frames; skipping.\n");
      return Action::SKIP;
    }
     mprintf("\tUsing existing temperature information in Frames.\n");
  } else {
    // Calculate temperature from velocities.
    if (cInfo_.HasTemp() && mode_ == CALC_AND_MODIFY)
      mprintf("Warning: Overwriting temperature information in Frames.\n");
    // Masks
    if (setup.Top().SetupIntegerMask( Mask_ )) return Action::ERR;
    Mask_.MaskInfo();
    if (Mask_.None()) {
      mprintf("Warning: temperature: No atoms selected in [%s]\n", Mask_.MaskString());
      return Action::SKIP;
    }
    // Calculate degrees of freedom taking into account constraints
    if (cons_.SetupConstraints( Mask_, setup.Top() )) return Action::ERR;
    // Determine d.o.f. offset if necessary
    dof_offset_ = 0;
    if (removeTrans_ && removeRot_) {
      if (Mask_.Nselected() == 1)
        dof_offset_ = 3;
      else if (Mask_.Nselected() == 2)
        dof_offset_ = 5;
      else
        dof_offset_ = 6;
    } else if (removeTrans_ || removeRot_) {
      dof_offset_ = 3;
    }
    if (dof_offset_ > 0) mprintf("\tRemoved %i additional degrees of freedom.\n", dof_offset_);
    // Update coordinate info if necessary.
    if (mode_ == CALC_AND_MODIFY) {
      cInfo_.SetTemperature( true );
      setup.SetCoordInfo( &cInfo_ );
      ret = Action::MODIFY_TOPOLOGY;
    }
  }
  return ret;
}

// Action_Temperature::DoAction()
Action::RetType Action_Temperature::DoAction(int frameNum, ActionFrame& frm) {
  double tdata;
  Action::RetType ret = Action::OK;
  if (mode_ == FROM_FRAME)
    tdata = frm.Frm().Temperature();
  else {
    tdata = frm.Frm().CalcTemperature(Mask_, cons_.DegreesOfFreedom() - dof_offset_);
    if (mode_ == CALC_AND_MODIFY) {
      frm.ModifyFrm().SetTemperature( tdata );
      ret = MODIFY_COORDS;
    }
  }
  Tdata_->Add(frameNum, &tdata);
  return ret;
}
