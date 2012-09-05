#include <cfloat>
#include "Action_Bounds.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Bounds::Action_Bounds() {}

/** Usage: bounds [<mask>] [out <filename>] */
int Action_Bounds::init() {
  outfilename_ = actionArgs.GetStringKey("out");
  
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  min_[0] = DBL_MAX;
  min_[1] = min_[0];
  min_[2] = min_[0];
  max_[0] = -DBL_MAX;
  max_[1] = max_[0];
  max_[2] = max_[0];

  mprintf("    BOUNDS: Calculating bounds for atoms in mask [%s]\n", mask_.MaskString());
  if (!outfilename_.empty())
    mprintf("\tOutput to file %s\n", outfilename_.c_str());
  return 0;
}

int Action_Bounds::setup() {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return 1;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: bounds: No atoms selected in mask.\n");
    return 1;
  }
  return 0;
}

int Action_Bounds::action() {
  for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
  {
    const double* xyz = currentFrame->XYZ( *atom );
    if (xyz[0] < min_[0]) min_[0] = xyz[0];
    if (xyz[0] > max_[0]) max_[0] = xyz[0];
    if (xyz[1] < min_[1]) min_[1] = xyz[1];
    if (xyz[1] > max_[1]) max_[1] = xyz[1];
    if (xyz[2] < min_[2]) min_[2] = xyz[2];
    if (xyz[2] > max_[2]) max_[2] = xyz[2];
  }
  return 0;
}

void Action_Bounds::print() {
  CpptrajFile outfile;

  if ( outfile.OpenWrite( outfilename_ ) ) return;
  outfile.Printf("%f < X < %f\n", min_[0], max_[0]);
  outfile.Printf("%f < Y < %f\n", min_[1], max_[1]);
  outfile.Printf("%f < Z < %f\n", min_[2], max_[2]);
  outfile.CloseFile();
}


