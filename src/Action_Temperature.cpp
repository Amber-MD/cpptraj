#include "Action_Temperature.h"
#include "CpptrajStdio.h"

Action_Temperature::Action_Temperature() :
  Tdata_(0),
  getTempFromFrame_(false),
  removeTrans_(false),
  removeRot_(false),
  dof_offset_(0)
{}

void Action_Temperature::Help() const {
  mprintf("\t[<name>] {frame | [<mask>] %s [remove {trans|rot|both}]}\n"
          "\t[out <filename>]\n", Constraints::constraintArgs);
  mprintf("  Calculate temperature in frame based on velocity information.\n"
          "  If 'frame' is specified just use frame temperature (read in from\n"
          "  e.g. REMD trajectory). The 'ntc' keyword can be used to correct\n"
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
    getTempFromFrame_ = true;
  else {
    getTempFromFrame_ = false;
    if (cons_.InitConstraints( actionArgs )) return Action::ERR;
  }
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Take into account removal of global degrees of freedom?
  removeTrans_ = false;
  removeRot_ = false;
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
  if (!getTempFromFrame_) {
    if (Mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  }
  // DataSet 
  Tdata_ =  init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "Tdata");
  if (Tdata_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddDataSet( Tdata_ );
  
  if (getTempFromFrame_) {
    mprintf("    TEMPERATURE: Frame temperatures will be saved in data set %s\n",
             Tdata_->legend());
  } else {
    mprintf("    TEMPERATURE: Calculate temperature for atoms in mask [%s]\n", Mask_.MaskString());
    mprintf("\tConstraints: %s\n", cons_.shakeString());
    if (removeTrans_) mprintf("\tAssuming translational degs. of freedom removed.\n");
    if (removeRot_) mprintf("\tAssuming rotational degs. of freedom removed.\n");
  }
  return Action::OK;
}

// Action_Temperature::Setup()
Action::RetType Action_Temperature::Setup(ActionSetup& setup) {
  if (getTempFromFrame_) {
    // Get temperature from frame.
    if (!setup.CoordInfo().HasTemperature()) {
      mprintf("Warning: No temperature information; skipping.\n");
      return Action::SKIP;
    }
  } else {
    // Calculate temperature from velocities.
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
  }
  return Action::OK;
}

// Action_Temperature::DoAction()
Action::RetType Action_Temperature::DoAction(int frameNum, ActionFrame& frm) {
  double tdata;
  if (getTempFromFrame_)
    tdata = frm.Frm().Temperature();
  else
    tdata = frm.Frm().CalcTemperature(Mask_, cons_.DegreesOfFreedom() - dof_offset_);
  Tdata_->Add(frameNum, &tdata);
  return Action::OK;
}
