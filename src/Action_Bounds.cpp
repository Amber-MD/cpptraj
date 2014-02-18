#include <cfloat>
#include "Action_Bounds.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Bounds::Action_Bounds() : ensembleNum_(-1) {}

void Action_Bounds::Help() {
  mprintf("\t[<mask>] [out <filename>] [dx <dx>] [dy <dy>] [dz <dz>]\n"
          "  Calcuate the max/min coordinates (X,Y,Z) of atoms in <mask>.\n");
}

Action::RetType Action_Bounds::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  outfilename_ = actionArgs.GetStringKey("out");
  dxyz_[0] = actionArgs.getKeyDouble("dx", -1.0);
  dxyz_[1] = actionArgs.getKeyDouble("dy", -1.0);
  dxyz_[2] = actionArgs.getKeyDouble("dz", -1.0);
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
  return Action::OK;
}

Action::RetType Action_Bounds::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: bounds: No atoms selected in mask.\n");
    return Action::ERR;
  }
  return Action::OK;
}

Action::RetType Action_Bounds::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
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
  return Action::OK;
}

static const char cXYZ[3] = {'X', 'Y', 'Z'};
void Action_Bounds::Print() {
  CpptrajFile outfile;
  if ( outfile.OpenEnsembleWrite( outfilename_, ensembleNum_ ) ) return;
  for (int i = 0; i < 3; i++) {
    outfile.Printf("%f < %c < %f", min_[i], cXYZ[i], max_[i]);
    if (dxyz_[i] > 0.0)
      outfile.Printf("\tCenter= %f  Bins=%f", (max_[i] + min_[i]) / 2.0,
                     (max_[i] - min_[i]) / dxyz_[i]);
    outfile.Printf("\n");
  }
  outfile.CloseFile();
}
