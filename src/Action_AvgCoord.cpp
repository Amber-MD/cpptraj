#include <cmath>
#include "Action_AvgCoord.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_AvgCoord::Action_AvgCoord() :
  calcMagnitude_(false),
  useMass_(false)
{ } 

void Action_AvgCoord::Help() {

}

// DESTRUCTOR
Action_AvgCoord::~Action_AvgCoord() {
  outfile_.CloseFile();
}

// Action_AvgCoord::init()
/** Expected call: avgcoord [<mask>] [mass] outfile <file> [magnitude]
  */
Action::RetType Action_AvgCoord::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  useMass_ = actionArgs.hasKey("mass");
  calcMagnitude_ = actionArgs.hasKey("magnitude");
  ArgList::ConstArg filename = actionArgs.getKeyString("outfile");
  if (filename==NULL) {
    mprinterr("Error: avgcoord: must specify output file with 'outfile'\n");
    return Action::ERR;
  }
  if (outfile_.SetupWrite( filename, debugIn )) return Action::ERR;
  if (outfile_.OpenFile()) return Action::ERR;

  // Get Masks
  Mask_.SetMaskString( actionArgs.getNextMask() );

  mprintf("    AVGCOORD: Average of atoms in mask [%s] will be output to %s\n",
          Mask_.MaskString(), filename);
  if (useMass_)
    mprintf("              Average will be mass-weighted.\n");
  if (calcMagnitude_)
    mprintf("              Magnitude of average will also be printed.\n");

  return Action::OK;
}

// Action_AvgCoord::setup()
Action::RetType Action_AvgCoord::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask(Mask_) ) return Action::ERR;
  if (Mask_.None()) {
    mprintf("Warning: AvgCoord::setup: Mask contains 0 atoms.\n");
    return Action::ERR;
  }
  mprintf("\t%s (%i atoms)\n",Mask_.MaskString(), Mask_.Nselected());

  return Action::OK;  
}

// Action_AvgCoord::action()
/** Calc avg of coordinates in frame. */
Action::RetType Action_AvgCoord::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double position[3];
 
  if (useMass_)
    currentFrame->CenterOfMass(position, Mask_);
  else 
    currentFrame->GeometricCenter(position, Mask_);
  // Calculate the magnitude
  if (calcMagnitude_) {
    double r2 = (position[0]*position[0]) + (position[1]*position[1]) + (position[2]*position[2]);
    double r = sqrt( r2 );
    outfile_.Printf("%i %lf %lf %lf %lf\n", frameNum+1, 
                        position[0], position[1], position[2], r);
  } else {
    outfile_.Printf("%i %lf %lf %lf\n",frameNum+1, position[0], position[1], position[2]);
  }

  return Action::OK;
} 

