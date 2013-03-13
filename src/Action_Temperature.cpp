#include "Action_Temperature.h"
#include "CpptrajStdio.h"

Action_Temperature::Action_Temperature() :
  Tdata_(0),
  ntc_(1)
{}

void Action_Temperature::Help() {
  mprintf("\t[<name>] [<mask>] [ntc <#>] [out <filename>]\n");
}

// Action_Temperature::Init()
Action::RetType Action_Temperature::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Keywords
  ntc_ = actionArgs.getNextInteger(1);
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Masks
  Mask_.SetMaskString( actionArgs.GetMaskNext() );
  // DataSet 
  Tdata_ =  DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "Tdata");
  if (Tdata_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddSet( Tdata_ );

  mprintf("    TEMPERATURE: Calculate temperature for atoms in mask [%s]\n", Mask_.MaskString());
  mprintf("\tUsing SHAKE (ntc) value of %i\n", ntc_);
  return Action::OK;
}

// Action_Temperature::Setup()
Action::RetType Action_Temperature::Setup(Topology* currentParm, Topology** parmAddress)
{
  // Masks
  if (currentParm->SetupIntegerMask( Mask_ )) return Action::ERR;
  Mask_.MaskInfo();
  if (Mask_.None()) {
    mprintf("Warning: temperature: No atoms selected in [%s]\n", Mask_.MaskString());
    return Action::ERR;
  }
  return Action::OK;
}

// Action_Temperature::DoAction()
Action::RetType Action_Temperature::DoAction(int frameNum, Frame* currentFrame, 
                                             Frame** frameAddress) 
{
  double tdata = currentFrame->Temperature(Mask_);
  Tdata_->Add(frameNum, &tdata);
  return Action::OK;
}
