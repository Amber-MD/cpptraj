#include "Action_Channel.h"
#include "CpptrajStdio.h"

Action_Channel::Action_Channel() {}

void Action_Channel::Help() {
  mprintf("channel");
  Grid::Help();
  mprintf("\n\t<solute mask> [<solvent mask>]\n");
}

Action::RetType Action_Channel::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                                     DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  // Get solute grid options
  if (solute_.GridInit( "CHANNEL", actionArgs ))
    return Action::ERR;
  // Make solvent grid a copy of solute
  solvent_ = solute_;
  // Get solute mask
  std::string maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty()) {
    mprinterr("Error: No solute mask specified.\n");
    return Action::ERR;
  }
  soluteMask_.SetMaskString( maskexpr );
  // Solvent mask
  maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty())
    maskexpr = ":WAT@O";
  solventMask_.SetMaskString( maskexpr );

  mprintf("    CHANNEL: Solute mask [%s], solvent mask [%s]\n",
          soluteMask_.MaskString(), solventMask_.MaskString());
  return Action::OK;
}

Action::RetType Action_Channel::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup grid, checks box info.
  if (solute_.GridSetup( currentParm )) return Action::ERR; 
  if (solvent_.GridSetup( currentParm )) return Action::ERR;
  // Setup masks
  if (currentParm->SetupIntegerMask( soluteMask_ ) ||
      currentParm->SetupIntegerMask( solventMask_)   ) 
    return Action::ERR;
  soluteMask_.MaskInfo();
  if (soluteMask_.None()) {
    mprinterr("Error: channel: No solute atoms selected.\n");
    return Action::ERR;
  }
  solventMask_.MaskInfo();
  if (solventMask_.None()) {
    mprinterr("Error: channel: No solvent atoms selected.\n");
    return Action::ERR;
  }
  return Action::OK;
}

Action::RetType Action_Channel::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) 
{
  solute_.GridFrame( *currentFrame, soluteMask_ );
  solvent_.GridFrame( *currentFrame, solventMask_ );
  return Action::OK;
}
