#include "Action_Temperature.h"
#include "CpptrajStdio.h"

Action_Temperature::Action_Temperature() :
  Tdata_(0),
  getTempFromFrame_(false)
{}

void Action_Temperature::Help() const {
  mprintf("\t[<name>] {frame | [<mask>] %s} [out <filename>]\n", Constraints::constraintArgs);
  mprintf("  Calculate temperature in frame based on velocity information.\n"
          "  If 'frame' is specified just use frame temperature (read in from\n"
          "  e.g. REMD trajectory)\n");
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
  // Masks
  if (!getTempFromFrame_)
    Mask_.SetMaskString( actionArgs.GetMaskNext() );
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
  }
  return Action::OK;
}

// Action_Temperature::Setup()
Action::RetType Action_Temperature::Setup(ActionSetup& setup) {
  if (!getTempFromFrame_) {
    // Masks
    if (setup.Top().SetupIntegerMask( Mask_ )) return Action::ERR;
    Mask_.MaskInfo();
    if (Mask_.None()) {
      mprintf("Warning: temperature: No atoms selected in [%s]\n", Mask_.MaskString());
      return Action::SKIP;
    }
    // Calculate degrees of freedom taking into account constraints
    if (cons_.SetupConstraints( Mask_, setup.Top() )) return Action::ERR;
  }
  return Action::OK;
}

// Action_Temperature::DoAction()
Action::RetType Action_Temperature::DoAction(int frameNum, ActionFrame& frm) {
  double tdata;
  if (getTempFromFrame_)
    tdata = frm.Frm().Temperature();
  else
    tdata = frm.Frm().CalcTemperature(Mask_, cons_.DegreesOfFreedom());
  Tdata_->Add(frameNum, &tdata);
  return Action::OK;
}
