#include <cmath>
#include "Action_AvgCoord.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_AvgCoord::Action_AvgCoord() :
  calcMagnitude_(false)
{
  useMass_ = false;
} 

// DESTRUCTOR
Action_AvgCoord::~Action_AvgCoord() {
  outfile_.CloseFile();
}

// Action_AvgCoord::init()
/** Expected call: avgcoord <mask> [mass] outfile <file> [magnitude]
  */
int Action_AvgCoord::init() {
  // Get keywords
  useMass_ = actionArgs.hasKey("mass");
  calcMagnitude_ = actionArgs.hasKey("magnitude");
  char *filename = actionArgs.getKeyString("outfile",NULL);
  if (filename==NULL) {
    mprinterr("Error: avgcoord: must specify output file with 'outfile'\n");
    return 1;
  }
  if (outfile_.SetupWrite( filename, debug )) return 1;
  if (outfile_.OpenFile()) return 1;

  // Get Masks
  char* mask1 = actionArgs.getNextMask();
  Mask_.SetMaskString(mask1);

  mprintf("    AVGCOORD: Average of atoms in mask [%s] will be output to %s\n",
          Mask_.MaskString(), filename);
  if (useMass_)
    mprintf("              Average will be mass-weighted.\n");
  if (calcMagnitude_)
    mprintf("              Magnitude of average will also be printed.\n");

  return 0;
}

// Action_AvgCoord::setup()
int Action_AvgCoord::setup() {
  if ( currentParm->SetupIntegerMask(Mask_) ) return 1;
  if (Mask_.None()) {
    mprintf("Warning: AvgCoord::setup: Mask contains 0 atoms.\n");
    return 1;
  }
  mprintf("\t%s (%i atoms)\n",Mask_.MaskString(), Mask_.Nselected());

  return 0;  
}

// Action_AvgCoord::action()
/** Calc avg of coordinates in frame. */
int Action_AvgCoord::action() {
  double position[3];
  
  currentFrame->GeometricCenter(&Mask_, position);
  // Calculate the magnitude
  if (calcMagnitude_) {
    double r2 = (position[0]*position[0]) + (position[1]*position[1]) + (position[2]*position[2]);
    double r = sqrt( r2 );
    outfile_.Printf("%i %lf %lf %lf %lf\n", frameNum+1, 
                        position[0], position[1], position[2], r);
  } else {
    outfile_.Printf("%i %lf %lf %lf\n",frameNum+1, position[0], position[1], position[2]);
  }

  return 0;
} 

